//! Database Query Optimizations for Genomic Data
//! =============================================
//!
//! Optimizes SQLite operations for high-throughput genomic workloads:
//! - Batch processing with prepared statements
//! - Optimal indexing strategies for genomic queries
//! - Connection pooling for concurrent access
//! - Query result caching and compression

use anyhow::{Context, Result};
use rusqlite::{params, Connection};
use serde_json;
use std::sync::Arc;
use parking_lot::{Mutex, RwLock};
use lz4_flex::{compress, decompress};
use dashmap::DashMap;
use serde::{Deserialize, Serialize};

/// High-performance database interface for genomic data
pub struct GenomicDatabase {
    /// Connection pool for concurrent access
    connection_pool: Arc<ConnectionPool>,
    /// Query result cache
    query_cache: Arc<DashMap<String, CachedResult>>,
    /// Configuration
    config: DatabaseConfig,
    /// Performance metrics
    metrics: Arc<DatabaseMetrics>,
}

#[derive(Debug, Clone)]
pub struct DatabaseConfig {
    /// Maximum connections in pool
    pub max_connections: usize,
    /// Batch size for bulk operations
    pub batch_size: usize,
    /// Enable query result caching
    pub enable_caching: bool,
    /// Cache size limit in MB
    pub cache_size_mb: usize,
    /// Enable compression for large results
    pub enable_compression: bool,
    /// WAL mode configuration
    pub wal_checkpoint_interval: usize,
}

impl Default for DatabaseConfig {
    fn default() -> Self {
        Self {
            max_connections: 8,
            batch_size: 10000,
            enable_caching: true,
            cache_size_mb: 512,
            enable_compression: true,
            wal_checkpoint_interval: 1000,
        }
    }
}

/// Connection pool for managing database connections
struct ConnectionPool {
    connections: Mutex<Vec<Connection>>,
    max_size: usize,
    database_path: String,
}

impl ConnectionPool {
    fn new(database_path: &str, max_size: usize) -> Result<Self> {
        let mut connections = Vec::with_capacity(max_size);
        
        for _ in 0..max_size {
            let mut conn = Connection::open(database_path)
                .context("Failed to create database connection")?;
            
            Self::configure_connection(&mut conn)?;
            connections.push(conn);
        }
        
        Ok(Self {
            connections: Mutex::new(connections),
            max_size,
            database_path: database_path.to_string(),
        })
    }
    
    fn configure_connection(conn: &mut Connection) -> Result<()> {
        // Handle all PRAGMA statements that might return results by using query pattern
        let pragma_statements = [
            "PRAGMA journal_mode = WAL",
            "PRAGMA cache_size = 20000", // ~80MB cache  
            "PRAGMA temp_store = memory",
            "PRAGMA mmap_size = 1073741824", // 1GB mmap
            "PRAGMA foreign_keys = ON", 
            "PRAGMA synchronous = NORMAL",
        ];
        
        for pragma in &pragma_statements {
            let mut stmt = conn.prepare(pragma)?;
            let mut rows = stmt.query_map([], |_| Ok(()))?;
            // Consume any result rows - some PRAGMA statements return values
            while let Some(_) = rows.next() {}
        }
        
        Ok(())
    }
    
    fn get_connection(&self) -> Result<Connection> {
        let mut connections = self.connections.lock();
        
        if let Some(conn) = connections.pop() {
            Ok(conn)
        } else {
            // Create new connection if pool is empty
            let mut conn = Connection::open(&self.database_path)?;
            Self::configure_connection(&mut conn)?;
            Ok(conn)
        }
    }
    
    fn return_connection(&self, conn: Connection) {
        let mut connections = self.connections.lock();
        if connections.len() < self.max_size {
            connections.push(conn);
        }
        // Otherwise, let it drop and close
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct CachedResult {
    data: Vec<u8>,
    compressed: bool,
    original_size: usize,
    created_at: std::time::SystemTime,
    hit_count: usize,
}

#[derive(Debug, Default)]
struct DatabaseMetrics {
    queries_executed: std::sync::atomic::AtomicUsize,
    cache_hits: std::sync::atomic::AtomicUsize,
    cache_misses: std::sync::atomic::AtomicUsize,
    batch_operations: std::sync::atomic::AtomicUsize,
    total_query_time_ms: std::sync::atomic::AtomicUsize,
}

impl GenomicDatabase {
    /// Create new database instance with optimization
    pub fn new(database_path: &str, config: DatabaseConfig) -> Result<Self> {
        let connection_pool = Arc::new(ConnectionPool::new(database_path, config.max_connections)?);
        
        // Initialize database schema if needed
        {
            let conn = connection_pool.get_connection()?;
            Self::initialize_schema(&conn)?;
            Self::create_optimized_indexes(&conn)?;
            connection_pool.return_connection(conn);
        }
        
        Ok(Self {
            connection_pool,
            query_cache: Arc::new(DashMap::new()),
            config,
            metrics: Arc::new(DatabaseMetrics::default()),
        })
    }
    
    /// Initialize database schema with optimizations
    fn initialize_schema(conn: &Connection) -> Result<()> {
        // K-mer table with optimal column types
        conn.execute(
            "CREATE TABLE IF NOT EXISTS kmers (
                id INTEGER PRIMARY KEY,
                sequence TEXT NOT NULL,
                hash BLOB NOT NULL,
                count INTEGER NOT NULL DEFAULT 1,
                first_seen TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )?;
        
        // Sequence table with compression support
        conn.execute(
            "CREATE TABLE IF NOT EXISTS sequences (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                sequence_data BLOB NOT NULL,
                length INTEGER NOT NULL,
                gc_content REAL,
                compressed BOOLEAN DEFAULT FALSE,
                source TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )?;
        
        // Graph edges table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS graph_edges (
                id INTEGER PRIMARY KEY,
                from_kmer_id INTEGER NOT NULL,
                to_kmer_id INTEGER NOT NULL,
                weight REAL NOT NULL DEFAULT 1.0,
                confidence REAL NOT NULL DEFAULT 1.0,
                edge_type TEXT DEFAULT 'overlap',
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (from_kmer_id) REFERENCES kmers(id),
                FOREIGN KEY (to_kmer_id) REFERENCES kmers(id)
            )",
            [],
        )?;
        
        // Assembly results table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS assembly_results (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL UNIQUE,
                total_contigs INTEGER NOT NULL,
                total_length INTEGER NOT NULL,
                n50 INTEGER,
                mean_coverage REAL,
                assembly_data BLOB,
                compressed BOOLEAN DEFAULT FALSE,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )?;
        
        Ok(())
    }
    
    /// Create optimized indexes for genomic queries
    fn create_optimized_indexes(conn: &Connection) -> Result<()> {
        let indexes = [
            // Primary lookup indexes
            "CREATE INDEX IF NOT EXISTS idx_kmers_sequence ON kmers(sequence)",
            "CREATE INDEX IF NOT EXISTS idx_kmers_hash ON kmers(hash)",
            "CREATE INDEX IF NOT EXISTS idx_kmers_count ON kmers(count DESC)",
            
            // Sequence lookup indexes
            "CREATE INDEX IF NOT EXISTS idx_sequences_name ON sequences(name)",
            "CREATE INDEX IF NOT EXISTS idx_sequences_length ON sequences(length)",
            "CREATE INDEX IF NOT EXISTS idx_sequences_gc ON sequences(gc_content)",
            
            // Graph traversal indexes
            "CREATE INDEX IF NOT EXISTS idx_edges_from ON graph_edges(from_kmer_id)",
            "CREATE INDEX IF NOT EXISTS idx_edges_to ON graph_edges(to_kmer_id)",
            "CREATE INDEX IF NOT EXISTS idx_edges_weight ON graph_edges(weight DESC)",
            
            // Composite indexes for common queries
            "CREATE INDEX IF NOT EXISTS idx_edges_from_to ON graph_edges(from_kmer_id, to_kmer_id)",
            "CREATE INDEX IF NOT EXISTS idx_kmers_count_seq ON kmers(count DESC, sequence)",
            
            // Time-based indexes
            "CREATE INDEX IF NOT EXISTS idx_kmers_last_updated ON kmers(last_updated)",
            "CREATE INDEX IF NOT EXISTS idx_sequences_created ON sequences(created_at)",
        ];
        
        for index_sql in &indexes {
            conn.execute(index_sql, [])?;
        }
        
        Ok(())
    }
    
    /// Batch insert k-mers with optimal performance
    pub fn batch_insert_kmers(&self, kmers: &[(String, Vec<u8>, u32)]) -> Result<usize> {
        let start_time = std::time::Instant::now();
        let conn = self.connection_pool.get_connection()?;
        
        let tx = conn.unchecked_transaction()?;
        
        // Prepare statement for batch insert
        let mut stmt = tx.prepare(
            "INSERT OR REPLACE INTO kmers (sequence, hash, count, last_updated) 
             VALUES (?1, ?2, ?3, CURRENT_TIMESTAMP)"
        )?;
        
        let mut inserted = 0;
        
        for chunk in kmers.chunks(self.config.batch_size) {
            for (sequence, hash, count) in chunk {
                stmt.execute(params![sequence, hash, count])?;
                inserted += 1;
            }
        }
        
        // Drop statement before committing transaction  
        drop(stmt);
        tx.commit()?;
        self.connection_pool.return_connection(conn);
        
        // Update metrics
        let elapsed_ms = start_time.elapsed().as_millis() as usize;
        self.metrics.batch_operations.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.metrics.total_query_time_ms.fetch_add(elapsed_ms, std::sync::atomic::Ordering::Relaxed);
        
        Ok(inserted)
    }
    
    /// Optimized k-mer frequency query with caching
    pub fn get_frequent_kmers(&self, min_count: u32, limit: Option<usize>) -> Result<Vec<KmerRecord>> {
        let cache_key = format!("frequent_kmers_{}_{}", min_count, limit.unwrap_or(0));
        
        // Check cache first
        if self.config.enable_caching {
            if let Some(mut cached) = self.query_cache.get_mut(&cache_key) {
                cached.hit_count += 1;
                self.metrics.cache_hits.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                
                let data = if cached.compressed {
                    decompress(&cached.data, cached.original_size)?
                } else {
                    cached.data.clone()
                };
                
                // For now, fall back to serde_json for simplicity
                return Ok(serde_json::from_slice(&data)?);
            }
        }
        
        // Execute query
        let conn = self.connection_pool.get_connection()?;
        let start_time = std::time::Instant::now();
        
        let query = if let Some(limit) = limit {
            format!(
                "SELECT id, sequence, hash, count FROM kmers 
                 WHERE count >= ? 
                 ORDER BY count DESC, sequence 
                 LIMIT {}",
                limit
            )
        } else {
            "SELECT id, sequence, hash, count FROM kmers 
             WHERE count >= ? 
             ORDER BY count DESC, sequence".to_string()
        };
        
        let mut stmt = conn.prepare(&query)?;
        let rows = stmt.query_map(params![min_count], |row| {
            Ok(KmerRecord {
                id: row.get(0)?,
                sequence: row.get(1)?,
                hash: row.get(2)?,
                count: row.get(3)?,
            })
        })?;
        
        let mut results = Vec::new();
        for row in rows {
            results.push(row?);
        }
        
        // Drop statement before returning connection to avoid borrowing conflicts
        drop(stmt);
        self.connection_pool.return_connection(conn);
        
        // Cache result
        if self.config.enable_caching && !results.is_empty() {
            // For now, use serde_json for simplicity
            let serialized = serde_json::to_vec(&results)?;
            let serialized_len = serialized.len();
            let (data, compressed) = if self.config.enable_compression && serialized_len > 1024 {
                (compress(&serialized), true)
            } else {
                (serialized, false)
            };
            
            self.query_cache.insert(cache_key, CachedResult {
                data,
                compressed,
                original_size: serialized_len,
                created_at: std::time::SystemTime::now(),
                hit_count: 0,
            });
        }
        
        // Update metrics
        let elapsed_ms = start_time.elapsed().as_millis() as usize;
        self.metrics.queries_executed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.metrics.cache_misses.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.metrics.total_query_time_ms.fetch_add(elapsed_ms, std::sync::atomic::Ordering::Relaxed);
        
        Ok(results)
    }
    
    /// Batch graph edge insertion with deduplication
    pub fn batch_insert_edges(&self, edges: &[(i64, i64, f64, f64)]) -> Result<usize> {
        let start_time = std::time::Instant::now();
        let conn = self.connection_pool.get_connection()?;
        
        let tx = conn.unchecked_transaction()?;
        
        // Use ON CONFLICT to handle duplicates efficiently
        let mut stmt = tx.prepare(
            "INSERT INTO graph_edges (from_kmer_id, to_kmer_id, weight, confidence) 
             VALUES (?1, ?2, ?3, ?4)
             ON CONFLICT(from_kmer_id, to_kmer_id) DO UPDATE SET
             weight = MAX(weight, ?3),
             confidence = MAX(confidence, ?4)"
        )?;
        
        let mut inserted = 0;
        
        for chunk in edges.chunks(self.config.batch_size) {
            for &(from_id, to_id, weight, confidence) in chunk {
                stmt.execute(params![from_id, to_id, weight, confidence])?;
                inserted += 1;
            }
        }
        
        drop(stmt);
        tx.commit()?;
        self.connection_pool.return_connection(conn);
        
        // Update metrics
        let elapsed_ms = start_time.elapsed().as_millis() as usize;
        self.metrics.batch_operations.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.metrics.total_query_time_ms.fetch_add(elapsed_ms, std::sync::atomic::Ordering::Relaxed);
        
        Ok(inserted)
    }
    
    /// Optimized graph traversal query
    pub fn get_graph_neighborhood(
        &self, 
        center_kmer_id: i64, 
        max_distance: usize
    ) -> Result<Vec<GraphEdge>> {
        let conn = self.connection_pool.get_connection()?;
        let start_time = std::time::Instant::now();
        
        // Use recursive CTE for efficient graph traversal
        let query = "
            WITH RECURSIVE graph_walk(node_id, distance, path) AS (
                -- Base case: start from the center node
                SELECT ?1 as node_id, 0 as distance, ?1 as path
                
                UNION ALL
                
                -- Recursive case: expand to neighbors
                SELECT 
                    ge.to_kmer_id,
                    gw.distance + 1,
                    gw.path || ',' || ge.to_kmer_id
                FROM graph_walk gw
                JOIN graph_edges ge ON gw.node_id = ge.from_kmer_id
                WHERE gw.distance < ?2
                  AND INSTR(gw.path, ',' || ge.to_kmer_id || ',') = 0  -- Avoid cycles
            )
            SELECT DISTINCT 
                ge.id, ge.from_kmer_id, ge.to_kmer_id, ge.weight, ge.confidence
            FROM graph_walk gw
            JOIN graph_edges ge ON (gw.node_id = ge.from_kmer_id OR gw.node_id = ge.to_kmer_id)
            WHERE gw.distance <= ?2
            ORDER BY ge.weight DESC";
        
        let results = {
            let mut stmt = conn.prepare(query)?;
            let rows = stmt.query_map(params![center_kmer_id, max_distance], |row| {
                Ok(GraphEdge {
                    id: row.get(0)?,
                    from_kmer_id: row.get(1)?,
                    to_kmer_id: row.get(2)?,
                    weight: row.get(3)?,
                    confidence: row.get(4)?,
                })
            })?;
            
            let mut results = Vec::new();
            for row in rows {
                results.push(row?);
            }
            results
        };
        
        self.connection_pool.return_connection(conn);
        
        // Update metrics
        let elapsed_ms = start_time.elapsed().as_millis() as usize;
        self.metrics.queries_executed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        self.metrics.total_query_time_ms.fetch_add(elapsed_ms, std::sync::atomic::Ordering::Relaxed);
        
        Ok(results)
    }
    
    /// Clear old cache entries to prevent memory bloat
    pub fn cleanup_cache(&self) {
        if !self.config.enable_caching {
            return;
        }
        
        let now = std::time::SystemTime::now();
        let cache_retention = std::time::Duration::from_secs(3600); // 1 hour
        
        self.query_cache.retain(|_, cached| {
            now.duration_since(cached.created_at).unwrap_or_default() < cache_retention
        });
        
        // Also limit cache size
        while self.estimate_cache_size_mb() > self.config.cache_size_mb {
            // Remove least recently used items
            if let Some(entry) = self.query_cache.iter().min_by_key(|entry| entry.value().hit_count) {
                let key = entry.key().clone();
                self.query_cache.remove(&key);
            } else {
                break;
            }
        }
    }
    
    /// Estimate cache size in MB
    fn estimate_cache_size_mb(&self) -> usize {
        self.query_cache.iter().map(|entry| entry.data.len()).sum::<usize>() / (1024 * 1024)
    }
    
    /// Get database performance metrics
    pub fn get_metrics(&self) -> DatabaseStats {
        DatabaseStats {
            queries_executed: self.metrics.queries_executed.load(std::sync::atomic::Ordering::Relaxed),
            cache_hits: self.metrics.cache_hits.load(std::sync::atomic::Ordering::Relaxed),
            cache_misses: self.metrics.cache_misses.load(std::sync::atomic::Ordering::Relaxed),
            batch_operations: self.metrics.batch_operations.load(std::sync::atomic::Ordering::Relaxed),
            total_query_time_ms: self.metrics.total_query_time_ms.load(std::sync::atomic::Ordering::Relaxed),
            cache_size_mb: self.estimate_cache_size_mb(),
            cache_entries: self.query_cache.len(),
        }
    }

    /// Insert a batch of KmerRecord objects (TDD-required method)
    /// This is a wrapper around batch_insert_kmers for compatibility with TDD tests
    pub fn insert_kmers_batch(&self, kmers: &[KmerRecord]) -> Result<usize> {
        // Convert KmerRecord to the format expected by batch_insert_kmers
        let converted_kmers: Vec<(String, Vec<u8>, u32)> = kmers
            .iter()
            .map(|k| (k.sequence.clone(), k.hash.clone(), k.count))
            .collect();
        
        self.batch_insert_kmers(&converted_kmers)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerRecord {
    pub id: i64,
    pub sequence: String,
    pub hash: Vec<u8>,
    pub count: u32,
}

#[derive(Debug, Clone)]
pub struct GraphEdge {
    pub id: i64,
    pub from_kmer_id: i64,
    pub to_kmer_id: i64,
    pub weight: f64,
    pub confidence: f64,
}

#[derive(Debug)]
pub struct DatabaseStats {
    pub queries_executed: usize,
    pub cache_hits: usize,
    pub cache_misses: usize,
    pub batch_operations: usize,
    pub total_query_time_ms: usize,
    pub cache_size_mb: usize,
    pub cache_entries: usize,
}

impl DatabaseStats {
    pub fn print_summary(&self) {
        println!("📊 Database Performance Statistics:");
        println!("   Queries executed: {}", self.queries_executed);
        println!("   Batch operations: {}", self.batch_operations);
        println!("   Total query time: {}ms", self.total_query_time_ms);
        
        if self.queries_executed > 0 {
            println!("   Average query time: {:.2}ms", 
                self.total_query_time_ms as f64 / self.queries_executed as f64);
        }
        
        let total_cache_requests = self.cache_hits + self.cache_misses;
        if total_cache_requests > 0 {
            let hit_rate = (self.cache_hits as f64 / total_cache_requests as f64) * 100.0;
            println!("   Cache hit rate: {:.1}%", hit_rate);
            println!("   Cache entries: {}", self.cache_entries);
            println!("   Cache size: {}MB", self.cache_size_mb);
        }
        
        if self.total_query_time_ms > 0 {
            let queries_per_second = (self.queries_executed * 1000) / self.total_query_time_ms;
            println!("   Throughput: {} queries/sec", queries_per_second);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    
    #[test]
    fn test_database_basic_operations() -> Result<()> {
        let temp_dir = tempdir()?;
        let db_path = temp_dir.path().join("test.db");
        
        let db = GenomicDatabase::new(
            db_path.to_str().unwrap(),
            DatabaseConfig::default(),
        )?;
        
        // Test batch k-mer insertion
        let kmers = vec![
            ("ATCG".to_string(), vec![1, 2, 3, 4], 5),
            ("GCTA".to_string(), vec![5, 6, 7, 8], 3),
            ("ATCG".to_string(), vec![1, 2, 3, 4], 8), // Duplicate with higher count
        ];
        
        let inserted = db.batch_insert_kmers(&kmers)?;
        assert_eq!(inserted, 3);
        
        // Test frequent k-mers query
        let frequent = db.get_frequent_kmers(3, Some(10))?;
        assert!(frequent.len() > 0);
        
        // Test edge insertion
        let edges = vec![
            (1, 2, 0.8, 0.9),
            (2, 3, 0.7, 0.8),
            (1, 3, 0.6, 0.7), // Should create triangle
        ];
        
        let edge_inserted = db.batch_insert_edges(&edges)?;
        assert_eq!(edge_inserted, 3);
        
        // Print statistics
        let stats = db.get_metrics();
        stats.print_summary();
        
        Ok(())
    }
    
    #[test]
    fn test_caching_performance() -> Result<()> {
        let temp_dir = tempdir()?;
        let db_path = temp_dir.path().join("cache_test.db");
        
        let db = GenomicDatabase::new(
            db_path.to_str().unwrap(),
            DatabaseConfig::default(),
        )?;
        
        // Insert test data
        let kmers = (0..1000).map(|i| {
            (format!("KMER{:04}", i), vec![i as u8; 4], (i % 100) as u32)
        }).collect::<Vec<_>>();
        
        db.batch_insert_kmers(&kmers)?;
        
        // First query (cache miss)
        let start = std::time::Instant::now();
        let result1 = db.get_frequent_kmers(50, Some(20))?;
        let first_query_time = start.elapsed();
        
        // Second query (cache hit)
        let start = std::time::Instant::now();
        let result2 = db.get_frequent_kmers(50, Some(20))?;
        let second_query_time = start.elapsed();
        
        // Results should be identical
        assert_eq!(result1.len(), result2.len());
        
        // Second query should be much faster
        println!("First query: {:?}", first_query_time);
        println!("Second query: {:?}", second_query_time);
        assert!(second_query_time < first_query_time / 2);
        
        let stats = db.get_metrics();
        assert!(stats.cache_hits > 0);
        stats.print_summary();
        
        Ok(())
    }
}