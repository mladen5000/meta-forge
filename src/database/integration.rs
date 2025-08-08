use std::path::{Path, PathBuf};
use ahash::{AHashMap, AHashSet};
use anyhow::{Result, Context};
use rusqlite::{Connection, params, OptionalExtension};
use serde::{Serialize, Deserialize};
use lz4_flex::{compress, decompress};
use dashmap::DashMap;
use parking_lot::RwLock;
use std::sync::Arc;

use crate::core::data_structures::{AssemblyStats, Contig};

/// Calculate GC content for a DNA sequence
fn calculate_gc_content(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }
    
    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C')
        .count();
    
    gc_count as f64 / sequence.len() as f64
}

/// Comprehensive database management for metagenomics pipeline
pub struct MetagenomicsDatabase {
    /// Primary database connection
    connection: Arc<RwLock<Connection>>,
    /// Database file path
    db_path: PathBuf,
    /// In-memory caches for performance
    caches: DatabaseCaches,
    /// Configuration
    config: DatabaseConfig,
}

#[derive(Debug, Clone)]
pub struct DatabaseConfig {
    /// Enable WAL mode for better concurrent access
    pub enable_wal_mode: bool,
    /// Cache size in pages
    pub cache_size: i32,
    /// Enable foreign keys
    pub enable_foreign_keys: bool,
    /// Batch size for bulk operations
    pub batch_size: usize,
    /// Enable compression for large data
    pub enable_compression: bool,
    /// Memory limit for caches (MB)
    pub cache_memory_limit_mb: usize,
}

impl Default for DatabaseConfig {
    fn default() -> Self {
        Self {
            enable_wal_mode: true,
            cache_size: 10000, // 10k pages (~40MB)
            enable_foreign_keys: true,
            batch_size: 1000,
            enable_compression: true,
            cache_memory_limit_mb: 256,
        }
    }
}

/// In-memory caches for frequently accessed data
struct DatabaseCaches {
    /// Taxonomy name cache
    taxonomy_cache: Arc<DashMap<u32, String>>,
    /// Sequence features cache
    features_cache: Arc<DashMap<String, CompressedFeatures>>,
    /// K-mer lookup cache
    kmer_cache: Arc<DashMap<String, Vec<SequenceMatch>>>,
    /// Assembly stats cache
    assembly_cache: Arc<DashMap<String, AssemblyRecord>>,
}

#[derive(Clone, Serialize, Deserialize)]
struct CompressedFeatures {
    compressed_data: Vec<u8>,
    feature_dim: usize,
    compression_ratio: f32,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SequenceMatch {
    pub sequence_id: i64,
    pub taxonomy_id: u32,
    pub similarity_score: f64,
    pub match_length: usize,
    pub e_value: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct AssemblyRecord {
    pub id: i64,
    pub sample_name: String,
    pub total_contigs: usize,
    pub total_length: usize,
    pub n50: usize,
    pub mean_coverage: f64,
    pub created_at: chrono::DateTime<chrono::Utc>,
    pub metadata: String,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct TaxonomyEntry {
    pub id: u32,
    pub name: String,
    pub lineage: String,
    pub rank: String,
    pub parent_id: Option<u32>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct SequenceEntry {
    pub id: i64,
    pub sequence_hash: String,
    pub sequence_data: String,
    pub length: usize,
    pub gc_content: f64,
    pub taxonomy_id: Option<u32>,
    pub source: String,
    pub created_at: chrono::DateTime<chrono::Utc>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct TrainingExample {
    pub id: i64,
    pub features: Vec<f64>,
    pub taxonomy_id: u32,
    pub sequence_id: i64,
    pub feature_version: String,
}

impl MetagenomicsDatabase {
    /// Create a new database instance
    pub fn new<P: AsRef<Path>>(db_path: P, config: DatabaseConfig) -> Result<Self> {
        let db_path = db_path.as_ref().to_path_buf();
        
        // Create directory if it doesn't exist
        if let Some(parent) = db_path.parent() {
            std::fs::create_dir_all(parent)
                .context("Failed to create database directory")?;
        }
        
        // Open database connection
        let conn = Connection::open(&db_path)
            .context("Failed to open database connection")?;
        
        let connection = Arc::new(RwLock::new(conn));
        
        let caches = DatabaseCaches {
            taxonomy_cache: Arc::new(DashMap::new()),
            features_cache: Arc::new(DashMap::new()),
            kmer_cache: Arc::new(DashMap::new()),
            assembly_cache: Arc::new(DashMap::new()),
        };
        
        let mut db = Self {
            connection,
            db_path,
            caches,
            config,
        };
        
        // Initialize database schema and settings
        db.initialize_database()?;
        
        println!("📁 Database initialized at {}", db.db_path.display());
        
        Ok(db)
    }
    
    /// Initialize database schema and settings
    fn initialize_database(&mut self) -> Result<()> {
        let conn = self.connection.write();
        
        // Configure database settings
        self.configure_database_settings(&conn)?;
        
        // Create all tables
        self.create_schema(&conn)?;
        
        // Create indices for performance
        self.create_indices(&conn)?;
        
        drop(conn);
        Ok(())
    }
    
    fn configure_database_settings(&self, conn: &Connection) -> Result<()> {
        // Enable WAL mode for better concurrent access
        if self.config.enable_wal_mode {
            conn.pragma_update(None, "journal_mode", "WAL")?;
        }
        
        // Set cache size
        conn.pragma_update(None, "cache_size", self.config.cache_size)?;
        
        // Enable foreign keys
        if self.config.enable_foreign_keys {
            conn.pragma_update(None, "foreign_keys", "ON")?;
        }
        
        // Other performance settings
        conn.pragma_update(None, "synchronous", "NORMAL")?;
        conn.pragma_update(None, "temp_store", "MEMORY")?;
        conn.pragma_update(None, "mmap_size", 268435456i64)?; // 256MB mmap
        
        Ok(())
    }
    
    fn create_schema(&self, conn: &Connection) -> Result<()> {
        // Taxonomy table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS taxonomy (
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL UNIQUE,
                lineage TEXT,
                rank TEXT,
                parent_id INTEGER,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (parent_id) REFERENCES taxonomy (id)
            )",
            [],
        )?;
        
        // Sequences table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS sequences (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence_hash TEXT NOT NULL UNIQUE,
                sequence_data TEXT,
                length INTEGER NOT NULL,
                gc_content REAL,
                taxonomy_id INTEGER,
                source TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (taxonomy_id) REFERENCES taxonomy (id)
            )",
            [],
        )?;
        
        // Features table (for ML training data)
        conn.execute(
            "CREATE TABLE IF NOT EXISTS sequence_features (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence_id INTEGER NOT NULL,
                features BLOB NOT NULL,
                feature_dim INTEGER NOT NULL,
                feature_version TEXT NOT NULL,
                taxonomy_id INTEGER,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (sequence_id) REFERENCES sequences (id),
                FOREIGN KEY (taxonomy_id) REFERENCES taxonomy (id)
            )",
            [],
        )?;
        
        // Assembly results table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS assemblies (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sample_name TEXT NOT NULL,
                total_contigs INTEGER NOT NULL,
                total_length INTEGER NOT NULL,
                longest_contig INTEGER,
                n50 INTEGER,
                n90 INTEGER,
                mean_coverage REAL,
                gc_content REAL,
                gaps INTEGER,
                metadata TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )?;
        
        // Contigs table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS contigs (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                assembly_id INTEGER NOT NULL,
                contig_name TEXT NOT NULL,
                sequence TEXT NOT NULL,
                length INTEGER NOT NULL,
                coverage REAL,
                gc_content REAL,
                annotations TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (assembly_id) REFERENCES assemblies (id)
            )",
            [],
        )?;
        
        // Annotations table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS annotations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                contig_id INTEGER,
                sequence_id INTEGER,
                taxonomy_id INTEGER,
                gene_name TEXT,
                gene_start INTEGER,
                gene_end INTEGER,
                strand TEXT,
                confidence REAL,
                method TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (contig_id) REFERENCES contigs (id),
                FOREIGN KEY (sequence_id) REFERENCES sequences (id),
                FOREIGN KEY (taxonomy_id) REFERENCES taxonomy (id)
            )",
            [],
        )?;
        
        // K-mer index table for fast lookups
        conn.execute(
            "CREATE TABLE IF NOT EXISTS kmer_index (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                kmer_hash TEXT NOT NULL,
                kmer_sequence TEXT NOT NULL,
                k_size INTEGER NOT NULL,
                sequence_id INTEGER NOT NULL,
                position INTEGER NOT NULL,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (sequence_id) REFERENCES sequences (id)
            )",
            [],
        )?;
        
        // Analysis results table
        conn.execute(
            "CREATE TABLE IF NOT EXISTS analysis_results (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                analysis_type TEXT NOT NULL,
                input_data TEXT NOT NULL,
                results BLOB NOT NULL,
                parameters TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            )",
            [],
        )?;
        
        println!("✅ Database schema created successfully");
        Ok(())
    }
    
    fn create_indices(&self, conn: &Connection) -> Result<()> {
        // Indices for fast lookups
        let indices = [
            "CREATE INDEX IF NOT EXISTS idx_sequences_hash ON sequences (sequence_hash);",
            "CREATE INDEX IF NOT EXISTS idx_sequences_taxonomy ON sequences (taxonomy_id);",
            "CREATE INDEX IF NOT EXISTS idx_sequences_length ON sequences (length);",
            "CREATE INDEX IF NOT EXISTS idx_features_sequence ON sequence_features (sequence_id);",
            "CREATE INDEX IF NOT EXISTS idx_features_taxonomy ON sequence_features (taxonomy_id);",
            "CREATE INDEX IF NOT EXISTS idx_features_version ON sequence_features (feature_version);",
            "CREATE INDEX IF NOT EXISTS idx_contigs_assembly ON contigs (assembly_id);",
            "CREATE INDEX IF NOT EXISTS idx_contigs_length ON contigs (length);",
            "CREATE INDEX IF NOT EXISTS idx_annotations_contig ON annotations (contig_id);",
            "CREATE INDEX IF NOT EXISTS idx_annotations_taxonomy ON annotations (taxonomy_id);",
            "CREATE INDEX IF NOT EXISTS idx_kmer_hash ON kmer_index (kmer_hash);",
            "CREATE INDEX IF NOT EXISTS idx_kmer_sequence ON kmer_index (sequence_id);",
            "CREATE INDEX IF NOT EXISTS idx_taxonomy_name ON taxonomy (name);",
            "CREATE INDEX IF NOT EXISTS idx_taxonomy_parent ON taxonomy (parent_id);",
            "CREATE INDEX IF NOT EXISTS idx_analysis_type ON analysis_results (analysis_type);",
        ];
        
        for index_sql in &indices {
            conn.execute(index_sql, [])?;
        }
        
        println!("✅ Database indices created successfully");
        Ok(())
    }
    
    /// Insert taxonomy entries
    pub fn insert_taxonomy_entries(&self, entries: &[TaxonomyEntry]) -> Result<()> {
        let conn = self.connection.write();
        let tx = conn.unchecked_transaction()?;
        
        {
            let mut stmt = tx.prepare(
                "INSERT OR IGNORE INTO taxonomy (id, name, lineage, rank, parent_id) 
                 VALUES (?1, ?2, ?3, ?4, ?5)"
            )?;
            
            for entry in entries {
                stmt.execute(params![
                    entry.id,
                    entry.name,
                    entry.lineage,
                    entry.rank,
                    entry.parent_id,
                ])?;
                
                // Update cache
                self.caches.taxonomy_cache.insert(entry.id, entry.name.clone());
            }
        } // stmt is dropped here
        
        tx.commit()?;
        println!("✅ Inserted {} taxonomy entries", entries.len());
        Ok(())
    }
    
    /// Insert sequence entries
    pub fn insert_sequences(&self, sequences: &[SequenceEntry]) -> Result<Vec<i64>> {
        let conn = self.connection.write();
        let tx = conn.unchecked_transaction()?;
        
        let mut sequence_ids = Vec::new();
        
        {
            let mut stmt = tx.prepare(
                "INSERT INTO sequences (sequence_hash, sequence_data, length, gc_content, taxonomy_id, source) 
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6)"
            )?;
            
            for seq in sequences {
                stmt.execute(params![
                    seq.sequence_hash,
                    seq.sequence_data,
                    seq.length,
                    seq.gc_content,
                    seq.taxonomy_id,
                    seq.source,
                ])?;
                
                sequence_ids.push(tx.last_insert_rowid());
            }
        } // stmt is dropped here
        
        tx.commit()?;
        println!("✅ Inserted {} sequences", sequences.len());
        Ok(sequence_ids)
    }
    
    /// Insert feature vectors for training
    pub fn insert_training_features(&self, examples: &[TrainingExample]) -> Result<()> {
        let conn = self.connection.write();
        let tx = conn.unchecked_transaction()?;
        
        {
            let mut stmt = tx.prepare(
                "INSERT INTO sequence_features (sequence_id, features, feature_dim, feature_version, taxonomy_id) 
                 VALUES (?1, ?2, ?3, ?4, ?5)"
            )?;
            
            for example in examples {
                // Compress features if enabled
                let features_data = if self.config.enable_compression {
                    let serialized = bincode::encode_to_vec(&example.features, bincode::config::standard())?;
                    compress(&serialized)
                } else {
                    bincode::encode_to_vec(&example.features, bincode::config::standard())?
                };
                
                stmt.execute(params![
                    example.sequence_id,
                    features_data,
                    example.features.len(),
                    example.feature_version,
                    example.taxonomy_id,
                ])?;
            }
        } // stmt is dropped here
        
        tx.commit()?;
        println!("✅ Inserted {} training examples", examples.len());
        Ok(())
    }
    
    /// Store assembly results
    pub fn store_assembly_results(&self, assembly: &AssemblyStats, sample_name: &str, metadata: &str) -> Result<i64> {
        let conn = self.connection.write();
        
        let mut stmt = conn.prepare(
            "INSERT INTO assemblies (sample_name, total_contigs, total_length, longest_contig, 
                                   n50, n90, mean_coverage, gc_content, gaps, metadata) 
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)"
        )?;
        
        stmt.execute(params![
            sample_name,
            assembly.num_contigs,
            assembly.total_length,
            assembly.largest_contig,
            assembly.n50,
            assembly.n90,
            assembly.coverage_mean,
            assembly.gc_content,
            0, // gaps field not available in AssemblyStats
            metadata,
        ])?;
        
        let assembly_id = conn.last_insert_rowid();
        
        // Cache the assembly record
        let record = AssemblyRecord {
            id: assembly_id,
            sample_name: sample_name.to_string(),
            total_contigs: assembly.num_contigs,
            total_length: assembly.total_length,
            n50: assembly.n50,
            mean_coverage: assembly.coverage_mean,
            created_at: chrono::Utc::now(),
            metadata: metadata.to_string(),
        };
        
        self.caches.assembly_cache.insert(sample_name.to_string(), record);
        
        println!("✅ Stored assembly results for sample: {sample_name}");
        Ok(assembly_id)
    }
    
    /// Store contigs from assembly
    pub fn store_contigs(&self, assembly_id: i64, contigs: &[Contig]) -> Result<()> {
        let conn = self.connection.write();
        let tx = conn.unchecked_transaction()?;
        
        {
            let mut stmt = tx.prepare(
                "INSERT INTO contigs (assembly_id, contig_name, sequence, length, coverage, gc_content) 
                 VALUES (?1, ?2, ?3, ?4, ?5, ?6)"
            )?;
            
            for contig in contigs {
                let contig_name = format!("contig_{}", contig.id);
                let gc_content = calculate_gc_content(&contig.sequence);
                
                stmt.execute(params![
                    assembly_id,
                    contig_name,
                    contig.sequence,
                    contig.length,
                    contig.coverage,
                    gc_content,
                ])?;
            }
        } // stmt is dropped here
        
        tx.commit()?;
        println!("✅ Stored {} contigs for assembly {}", contigs.len(), assembly_id);
        Ok(())
    }
    
    /// Query functions
    
    /// Get taxonomy name by ID (with caching)
    pub fn get_taxonomy_name(&self, taxonomy_id: u32) -> Result<Option<String>> {
        // Check cache first
        if let Some(name) = self.caches.taxonomy_cache.get(&taxonomy_id) {
            return Ok(Some(name.clone()));
        }
        
        // Query database
        let conn = self.connection.read();
        let mut stmt = conn.prepare("SELECT name FROM taxonomy WHERE id = ?1")?;
        
        let name: Option<String> = stmt.query_row([taxonomy_id], |row| {
            row.get(0)
        }).optional()?;
        
        // Cache result
        if let Some(ref name_str) = name {
            self.caches.taxonomy_cache.insert(taxonomy_id, name_str.clone());
        }
        
        Ok(name)
    }
    
    /// Get training data for machine learning
    pub fn get_training_data(&self, feature_version: &str, limit: Option<usize>) -> Result<Vec<TrainingExample>> {
        let conn = self.connection.read();
        
        let sql = if let Some(limit_val) = limit {
            format!(
                "SELECT sf.id, sf.features, sf.feature_dim, sf.taxonomy_id, sf.sequence_id 
                 FROM sequence_features sf 
                 WHERE sf.feature_version = ?1 
                 ORDER BY RANDOM() 
                 LIMIT {limit_val}"
            )
        } else {
            "SELECT sf.id, sf.features, sf.feature_dim, sf.taxonomy_id, sf.sequence_id 
             FROM sequence_features sf 
             WHERE sf.feature_version = ?1".to_string()
        };
        
        let mut stmt = conn.prepare(&sql)?;
        let rows = stmt.query_map([feature_version], |row| {
            let id: i64 = row.get(0)?;
            let features_blob: Vec<u8> = row.get(1)?;
            let feature_dim: usize = row.get(2)?;
            let taxonomy_id: u32 = row.get(3)?;
            let sequence_id: i64 = row.get(4)?;
            
            // Decompress features if necessary
            let features: Vec<f64> = if self.config.enable_compression {
                let decompressed = decompress(&features_blob, feature_dim * 8) // f64 = 8 bytes
                    .map_err(|e| rusqlite::Error::InvalidColumnType(1, "features".to_string(), rusqlite::types::Type::Blob))?;
                bincode::decode_from_slice(&decompressed, bincode::config::standard()).map(|(result, _)| result)
                    .map_err(|e| rusqlite::Error::InvalidColumnType(1, "features".to_string(), rusqlite::types::Type::Blob))?
            } else {
                bincode::decode_from_slice(&features_blob, bincode::config::standard()).map(|(result, _)| result)
                    .map_err(|e| rusqlite::Error::InvalidColumnType(1, "features".to_string(), rusqlite::types::Type::Blob))?
            };
            
            Ok(TrainingExample {
                id,
                features,
                taxonomy_id,
                sequence_id,
                feature_version: feature_version.to_string(),
            })
        })?;
        
        let mut examples = Vec::new();
        for row_result in rows {
            examples.push(row_result?);
        }
        
        println!("📊 Retrieved {} training examples", examples.len());
        Ok(examples)
    }
    
    /// Find similar sequences using k-mer matching
    pub fn find_similar_sequences(&self, query_sequence: &str, k: usize, max_results: usize) -> Result<Vec<SequenceMatch>> {
        // Generate cache key
        let cache_key = format!("{query_sequence}_{k}");
        
        // Check cache
        if let Some(cached_results) = self.caches.kmer_cache.get(&cache_key) {
            return Ok(cached_results.clone());
        }
        
        // Extract k-mers from query
        let query_kmers: AHashSet<String> = query_sequence
            .as_bytes()
            .windows(k)
            .filter_map(|window| std::str::from_utf8(window).ok())
            .map(|s| s.to_string())
            .collect();
        
        if query_kmers.is_empty() {
            return Ok(Vec::new());
        }
        
        let conn = self.connection.read();
        
        // Find sequences with matching k-mers
        let mut sequence_scores = AHashMap::new();
        
        for kmer in &query_kmers {
            let mut stmt = conn.prepare(
                "SELECT ki.sequence_id, s.taxonomy_id, s.length 
                 FROM kmer_index ki 
                 JOIN sequences s ON ki.sequence_id = s.id 
                 WHERE ki.kmer_sequence = ?1"
            )?;
            
            let rows = stmt.query_map([kmer], |row| {
                let sequence_id: i64 = row.get(0)?;
                let taxonomy_id: Option<u32> = row.get(1)?;
                let length: usize = row.get(2)?;
                Ok((sequence_id, taxonomy_id.unwrap_or(0), length))
            })?;
            
            for row_result in rows {
                let (sequence_id, taxonomy_id, length) = row_result?;
                let entry = sequence_scores.entry(sequence_id).or_insert((0, taxonomy_id, length));
                entry.0 += 1; // Increment match count
            }
        }
        
        // Calculate similarity scores and create results
        let mut results = Vec::new();
        for (sequence_id, (match_count, taxonomy_id, length)) in sequence_scores {
            let similarity_score = match_count as f64 / query_kmers.len() as f64;
            let e_value = self.calculate_e_value(match_count, query_kmers.len(), length);
            
            results.push(SequenceMatch {
                sequence_id,
                taxonomy_id,
                similarity_score,
                match_length: match_count,
                e_value,
            });
        }
        
        // Sort by similarity score (descending)
        results.sort_by(|a, b| b.similarity_score.partial_cmp(&a.similarity_score).unwrap());
        results.truncate(max_results);
        
        // Cache results
        self.caches.kmer_cache.insert(cache_key, results.clone());
        
        Ok(results)
    }
    
    /// Get assembly statistics
    pub fn get_assembly_stats(&self, sample_name: &str) -> Result<Option<AssemblyRecord>> {
        // Check cache first
        if let Some(cached) = self.caches.assembly_cache.get(sample_name) {
            return Ok(Some(cached.clone()));
        }
        
        let conn = self.connection.read();
        let mut stmt = conn.prepare(
            "SELECT id, sample_name, total_contigs, total_length, n50, mean_coverage, created_at, metadata 
             FROM assemblies WHERE sample_name = ?1"
        )?;
        
        let record = stmt.query_row([sample_name], |row| {
            Ok(AssemblyRecord {
                id: row.get(0)?,
                sample_name: row.get(1)?,
                total_contigs: row.get(2)?,
                total_length: row.get(3)?,
                n50: row.get(4)?,
                mean_coverage: row.get(5)?,
                created_at: row.get(6)?,
                metadata: row.get(7)?,
            })
        }).optional()?;
        
        // Cache the result
        if let Some(ref rec) = record {
            self.caches.assembly_cache.insert(sample_name.to_string(), rec.clone());
        }
        
        Ok(record)
    }
    
    /// Bulk operations for better performance
    
    /// Build k-mer index for fast sequence lookups
    pub fn build_kmer_index(&self, k_size: usize) -> Result<()> {
        println!("🔨 Building k-mer index with k={k_size}");
        
        let conn = self.connection.write();
        let tx = conn.unchecked_transaction()?;
        
        // Clear existing k-mer index for this k-size
        tx.execute("DELETE FROM kmer_index WHERE k_size = ?1", [k_size])?;
        
        // Collect sequences first
        let sequences = {
            let mut seq_stmt = tx.prepare("SELECT id, sequence_data FROM sequences WHERE length >= ?1")?;
            let sequences = seq_stmt.query_map([k_size], |row| {
                let id: i64 = row.get(0)?;
                let sequence: String = row.get(1)?;
                Ok((id, sequence))
            })?;
            sequences.collect::<Result<Vec<_>, _>>()?
        };
        
        // Insert k-mers
        let kmer_count = {
            let mut kmer_stmt = tx.prepare(
                "INSERT INTO kmer_index (kmer_hash, kmer_sequence, k_size, sequence_id, position) 
                 VALUES (?1, ?2, ?3, ?4, ?5)"
            )?;
            
            let mut count = 0;
            for (sequence_id, sequence) in sequences {
                for (pos, window) in sequence.as_bytes().windows(k_size).enumerate() {
                    if let Ok(kmer) = std::str::from_utf8(window) {
                        let kmer_hash = ahash::RandomState::new().hash_one(kmer);
                        
                        kmer_stmt.execute(params![
                            format!("{:x}", kmer_hash),
                            kmer,
                            k_size,
                            sequence_id,
                            pos,
                        ])?;
                        
                        count += 1;
                        
                        if count % 10000 == 0 {
                            println!("  Processed {count} k-mers");
                        }
                    }
                }
            }
            count
        }; // kmer_stmt is dropped here, return count
        
        tx.commit()?;
        println!("✅ Built k-mer index with {kmer_count} k-mers");
        Ok(())
    }
    
    /// Database maintenance operations
    
    /// Vacuum database to reclaim space
    pub fn vacuum_database(&self) -> Result<()> {
        println!("🧹 Vacuuming database...");
        let conn = self.connection.write();
        conn.execute("VACUUM;", [])?;
        println!("✅ Database vacuumed successfully");
        Ok(())
    }
    
    /// Analyze database for query optimization
    pub fn analyze_database(&self) -> Result<()> {
        println!("📊 Analyzing database for optimization...");
        let conn = self.connection.write();
        conn.execute("ANALYZE;", [])?;
        println!("✅ Database analysis completed");
        Ok(())
    }
    
    /// Get database statistics
    pub fn get_database_stats(&self) -> Result<DatabaseStats> {
        let conn = self.connection.read();
        
        let mut stats = DatabaseStats::default();
        
        // Count records in each table
        stats.taxonomy_count = conn.query_row("SELECT COUNT(*) FROM taxonomy", [], |row| row.get(0))?;
        stats.sequences_count = conn.query_row("SELECT COUNT(*) FROM sequences", [], |row| row.get(0))?;
        stats.features_count = conn.query_row("SELECT COUNT(*) FROM sequence_features", [], |row| row.get(0))?;
        stats.assemblies_count = conn.query_row("SELECT COUNT(*) FROM assemblies", [], |row| row.get(0))?;
        stats.contigs_count = conn.query_row("SELECT COUNT(*) FROM contigs", [], |row| row.get(0))?;
        stats.annotations_count = conn.query_row("SELECT COUNT(*) FROM annotations", [], |row| row.get(0))?;
        stats.kmer_index_count = conn.query_row("SELECT COUNT(*) FROM kmer_index", [], |row| row.get(0))?;
        
        // Database file size
        if let Ok(metadata) = std::fs::metadata(&self.db_path) {
            stats.db_size_bytes = metadata.len();
        }
        
        // Cache statistics
        stats.cache_hits = self.caches.taxonomy_cache.len() + 
                          self.caches.features_cache.len() + 
                          self.caches.kmer_cache.len() + 
                          self.caches.assembly_cache.len();
        
        Ok(stats)
    }
    
    /// Clear all caches
    pub fn clear_caches(&self) {
        self.caches.taxonomy_cache.clear();
        self.caches.features_cache.clear();
        self.caches.kmer_cache.clear();
        self.caches.assembly_cache.clear();
        println!("🗑️  All caches cleared");
    }
    
    // Helper methods
    
    fn calculate_e_value(&self, matches: usize, query_length: usize, subject_length: usize) -> f64 {
        // Simplified E-value calculation
        let m = query_length as f64;
        let n = subject_length as f64;
        let s = matches as f64;
        
        // Rough approximation
        (m * n) / (2.0_f64.powf(s))
    }
}

#[derive(Debug, Default)]
pub struct DatabaseStats {
    pub taxonomy_count: i64,
    pub sequences_count: i64,
    pub features_count: i64,
    pub assemblies_count: i64,
    pub contigs_count: i64,
    pub annotations_count: i64,
    pub kmer_index_count: i64,
    pub db_size_bytes: u64,
    pub cache_hits: usize,
}

impl std::fmt::Display for DatabaseStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "📊 Database Statistics:")?;
        writeln!(f, "  Taxonomy entries: {}", self.taxonomy_count)?;
        writeln!(f, "  Sequences: {}", self.sequences_count)?;
        writeln!(f, "  Feature vectors: {}", self.features_count)?;
        writeln!(f, "  Assemblies: {}", self.assemblies_count)?;
        writeln!(f, "  Contigs: {}", self.contigs_count)?;
        writeln!(f, "  Annotations: {}", self.annotations_count)?;
        writeln!(f, "  K-mer index entries: {}", self.kmer_index_count)?;
        writeln!(f, "  Database size: {:.2} MB", self.db_size_bytes as f64 / (1024.0 * 1024.0))?;
        writeln!(f, "  Cache entries: {}", self.cache_hits)?;
        Ok(())
    }
}

/// Database migration utilities
pub struct DatabaseMigrator {
    db: MetagenomicsDatabase,
}

impl DatabaseMigrator {
    pub fn new(db: MetagenomicsDatabase) -> Self {
        Self { db }
    }
    
    /// Import taxonomy data from standard files
    pub fn import_taxonomy_from_ncbi<P: AsRef<Path>>(&self, taxonomy_file: P) -> Result<()> {
        println!("📥 Importing NCBI taxonomy data...");
        
        let file_content = std::fs::read_to_string(taxonomy_file)?;
        let mut entries = Vec::new();
        
        for (line_num, line) in file_content.lines().enumerate() {
            if line_num == 0 { continue; } // Skip header
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 4 {
                let entry = TaxonomyEntry {
                    id: fields[0].parse()?,
                    name: fields[1].to_string(),
                    lineage: fields.get(2).unwrap_or(&"").to_string(),
                    rank: fields.get(3).unwrap_or(&"").to_string(),
                    parent_id: fields.get(4).and_then(|s| s.parse().ok()),
                };
                entries.push(entry);
            }
        }
        
        self.db.insert_taxonomy_entries(&entries)?;
        println!("✅ Imported {} taxonomy entries", entries.len());
        Ok(())
    }
    
    /// Import sequence data from FASTA files
    pub fn import_sequences_from_fasta<P: AsRef<Path> + std::fmt::Debug>(&self, fasta_file: P, source: &str) -> Result<()> {
        use bio::io::fasta;
        
        println!("📥 Importing sequences from FASTA...");
        
        let reader = fasta::Reader::from_file(fasta_file)?;
        let mut sequences = Vec::new();
        
        for (i, record_result) in reader.records().enumerate() {
            let record = record_result?;
            let sequence_str = std::str::from_utf8(record.seq())?;
            
            let sequence_hash = format!("{:x}", ahash::RandomState::new().hash_one(sequence_str));
            let gc_content = calculate_gc_content(sequence_str);
            
            let entry = SequenceEntry {
                id: 0, // Will be assigned by database
                sequence_hash,
                sequence_data: sequence_str.to_string(),
                length: sequence_str.len(),
                gc_content,
                taxonomy_id: None, // Would parse from header in real implementation
                source: source.to_string(),
                created_at: chrono::Utc::now(),
            };
            
            sequences.push(entry);
            
            if sequences.len() >= 1000 {
                self.db.insert_sequences(&sequences)?;
                sequences.clear();
            }
            
            if (i + 1) % 10000 == 0 {
                println!("  Processed {} sequences", i + 1);
            }
        }
        
        // Insert remaining sequences
        if !sequences.is_empty() {
            self.db.insert_sequences(&sequences)?;
        }
        
        println!("✅ Sequence import completed");
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    
    #[test]
    fn test_database_creation() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        
        let config = DatabaseConfig::default();
        let db = MetagenomicsDatabase::new(db_path, config).unwrap();
        
        let stats = db.get_database_stats().unwrap();
        assert_eq!(stats.taxonomy_count, 0);
        assert_eq!(stats.sequences_count, 0);
    }
    
    #[test]
    fn test_taxonomy_operations() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        
        let config = DatabaseConfig::default();
        let db = MetagenomicsDatabase::new(db_path, config).unwrap();
        
        let entries = vec![
            TaxonomyEntry {
                id: 1,
                name: "Escherichia coli".to_string(),
                lineage: "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia".to_string(),
                rank: "species".to_string(),
                parent_id: Some(562),
            },
        ];
        
        db.insert_taxonomy_entries(&entries).unwrap();
        
        let name = db.get_taxonomy_name(1).unwrap();
        assert_eq!(name, Some("Escherichia coli".to_string()));
        
        let stats = db.get_database_stats().unwrap();
        assert_eq!(stats.taxonomy_count, 1);
    }
    
    #[test]
    fn test_sequence_operations() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        
        let config = DatabaseConfig::default();
        let db = MetagenomicsDatabase::new(db_path, config).unwrap();
        
        // Insert required taxonomy entry first
        let taxonomy_entries = vec![
            TaxonomyEntry {
                id: 1,
                name: "Test Organism".to_string(),
                rank: "species".to_string(),
                parent_id: None,
                lineage: "cellular organisms; Bacteria; Test Organism".to_string(),
            },
        ];
        db.insert_taxonomy_entries(&taxonomy_entries).unwrap();
        
        let sequences = vec![
            SequenceEntry {
                id: 0,
                sequence_hash: "test_hash".to_string(),
                sequence_data: "ATCGATCGATCG".to_string(),
                length: 12,
                gc_content: 0.5,
                taxonomy_id: Some(1),
                source: "test".to_string(),
                created_at: chrono::Utc::now(),
            },
        ];
        
        let ids = db.insert_sequences(&sequences).unwrap();
        assert_eq!(ids.len(), 1);
        
        let stats = db.get_database_stats().unwrap();
        assert_eq!(stats.sequences_count, 1);
    }
    
    #[test]
    fn test_feature_compression() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");
        
        let mut config = DatabaseConfig::default();
        config.enable_compression = true;
        
        let db = MetagenomicsDatabase::new(db_path, config).unwrap();
        
        let features = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let examples = vec![
            TrainingExample {
                id: 0,
                features,
                taxonomy_id: 1,
                sequence_id: seq_ids[0], // Use actual sequence ID from database
                feature_version: "v1.0".to_string(),
            },
        ];
        
        db.insert_training_features(&examples).unwrap();
        
        let retrieved = db.get_training_data("v1.0", Some(10)).unwrap();
        assert_eq!(retrieved.len(), 1);
        assert_eq!(retrieved[0].features, examples[0].features);
    }
}