use anyhow::{anyhow, Result};
use memmap2::{MmapMut, MmapOptions};
use std::fs::{File, OpenOptions};
use std::mem::size_of;
use std::path::Path;
use std::ptr;
use std::sync::Arc; // kept if you need Arc elsewhere

use crate::core::data_structures::*; // AssemblyGraph, GraphNode, GraphEdge, CanonicalKmer, etc.

// ---------------------------- Public Views (safe API) ----------------------------

/// Public, copyable view of a node (safe to expose)
#[derive(Debug, Copy, Clone)]
pub struct NodeView {
    pub hash: u64,
    pub coverage: u32,
    pub degree: u32,
}

/// Public, copyable view of an edge (safe to expose)
#[derive(Debug, Copy, Clone)]
pub struct EdgeView {
    pub from_hash: u64,
    pub to_hash: u64,
    pub weight: u32,
    pub confidence: f32,
}

/// Memory usage statistics (public)
#[derive(Debug, Clone)]
pub struct MemoryStats {
    pub nodes: usize,
    pub edges: usize,
    pub node_memory_mb: u64,
    pub edge_memory_mb: u64,
    pub index_memory_mb: u64,
    pub total_memory_mb: u64,
}

impl MemoryStats {
    pub fn print_summary(&self) {
        println!("📊 Memory-Mapped Graph Statistics:");
        println!("   Nodes: {}", self.nodes);
        println!("   Edges: {}", self.edges);
        println!("   Node storage: {} MB", self.node_memory_mb);
        println!("   Edge storage: {} MB", self.edge_memory_mb);
        println!("   Index storage: {} MB", self.index_memory_mb);
        println!("   Total storage: {} MB", self.total_memory_mb);

        // Rough RAM savings estimate vs. a naive in-RAM structure (tunable)
        let estimated_ram_usage_mb =
            ((self.nodes as u64 * 64) + (self.edges as u64 * 32)) / (1024 * 1024);
        if estimated_ram_usage_mb > 0 {
            let savings_percent = if estimated_ram_usage_mb > self.total_memory_mb {
                (estimated_ram_usage_mb - self.total_memory_mb) as f64
                    / estimated_ram_usage_mb as f64
                    * 100.0
            } else {
                0.0
            };
            println!("   Estimated RAM savings: {:.1}%", savings_percent);
        }
    }
}

// ---------------------------- Internal POD (private) ----------------------------

// NOTE: these are kept private and never borrowed as &T from the mmap.
// We read/write them with *unaligned* ops, then convert to public views.

#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct MappedNode {
    hash: u64,
    coverage: u32,
    degree: u32,
    _padding: [u8; 4], // keep struct size a multiple of 8; size = 24
}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct MappedEdge {
    from_hash: u64,
    to_hash: u64,
    weight: u32,
    confidence: f32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug)]
struct IndexEntry {
    hash: u64,
    offset: u64,
}

// ------------------------------ Mapped Graph Core -------------------------------

pub struct MappedGraph {
    node_mmap: MmapMut,
    edge_mmap: MmapMut,
    index_mmap: MmapMut,
    node_count: usize,
    edge_count: usize,
    file_paths: MappedFilePaths,
}

struct MappedFilePaths {
    nodes_file: String,
    edges_file: String,
    index_file: String,
}

impl MappedGraph {
    /// Create a new memory-mapped graph with specified capacity
    pub fn new<P: AsRef<Path>>(
        base_path: P,
        estimated_nodes: usize,
        estimated_edges: usize,
    ) -> Result<Self> {
        let base = base_path.as_ref();

        // Sizes of internal structs
        let node_size = size_of::<MappedNode>();
        let edge_size = size_of::<MappedEdge>();
        let index_entry_size = size_of::<IndexEntry>();

        // File sizes based on capacities
        let nodes_file_size = estimated_nodes
            .checked_mul(node_size)
            .ok_or_else(|| anyhow!("nodes_file_size overflow"))?;
        let edges_file_size = estimated_edges
            .checked_mul(edge_size)
            .ok_or_else(|| anyhow!("edges_file_size overflow"))?;
        let index_file_size = estimated_nodes
            .checked_mul(index_entry_size)
            .ok_or_else(|| anyhow!("index_file_size overflow"))?;

        let file_paths = MappedFilePaths {
            nodes_file: base.join("nodes.mmap").to_string_lossy().to_string(),
            edges_file: base.join("edges.mmap").to_string_lossy().to_string(),
            index_file: base.join("index.mmap").to_string_lossy().to_string(),
        };

        // Ensure directory exists
        if let Some(parent) = base.parent() {
            std::fs::create_dir_all(parent)?;
        }

        // Create + size files
        let nodes_file = Self::create_sized_file(&file_paths.nodes_file, nodes_file_size)?;
        let edges_file = Self::create_sized_file(&file_paths.edges_file, edges_file_size)?;
        let index_file = Self::create_sized_file(&file_paths.index_file, index_file_size)?;

        // Map them
        let node_mmap = unsafe { MmapOptions::new().map_mut(&nodes_file)? };
        let edge_mmap = unsafe { MmapOptions::new().map_mut(&edges_file)? };
        let index_mmap = unsafe { MmapOptions::new().map_mut(&index_file)? };

        Ok(Self {
            node_mmap,
            edge_mmap,
            index_mmap,
            node_count: 0,
            edge_count: 0,
            file_paths,
        })
    }

    /// Add a node to the memory-mapped storage (unaligned write)
    pub fn add_node(&mut self, hash: u64, coverage: u32) -> Result<()> {
        let sz = size_of::<MappedNode>();
        if self.node_count * sz >= self.node_mmap.len() {
            return Err(anyhow!("Node storage full"));
        }
        let offset = self.node_count * sz;
        let node = MappedNode {
            hash,
            coverage,
            degree: 0,
            _padding: [0; 4],
        };
        unsafe {
            let p = self.node_mmap.as_mut_ptr().add(offset) as *mut MappedNode;
            ptr::write_unaligned(p, node);
        }
        self.update_index(hash, offset as u64)?;
        self.node_count += 1;
        Ok(())
    }

    /// Add an edge to the memory-mapped storage (unaligned write)
    pub fn add_edge(&mut self, from_hash: u64, to_hash: u64, weight: u32) -> Result<()> {
        let sz = size_of::<MappedEdge>();
        if self.edge_count * sz >= self.edge_mmap.len() {
            return Err(anyhow!("Edge storage full"));
        }
        let offset = self.edge_count * sz;
        let edge = MappedEdge {
            from_hash,
            to_hash,
            weight,
            confidence: 1.0,
        };
        unsafe {
            let p = self.edge_mmap.as_mut_ptr().add(offset) as *mut MappedEdge;
            ptr::write_unaligned(p, edge);
        }
        self.edge_count += 1;
        Ok(())
    }

    /// Get a node by hash (owned copy; no private type leak)
    pub fn get_node(&self, hash: u64) -> Option<NodeView> {
        self.lookup_node_offset(hash).map(|off| unsafe {
            let p = self.node_mmap.as_ptr().add(off as usize) as *const MappedNode;
            let n = ptr::read_unaligned(p);
            NodeView {
                hash: n.hash,
                coverage: n.coverage,
                degree: n.degree,
            }
        })
    }

    /// Iterate all nodes (opaque iterator; yields NodeView by value)
    pub fn iter_nodes(&self) -> impl Iterator<Item = NodeView> + '_ {
        let sz = size_of::<MappedNode>();
        let total = self.node_count;
        let mmap = &self.node_mmap;
        (0..total).map(move |i| unsafe {
            let off = i * sz;
            let p = mmap.as_ptr().add(off) as *const MappedNode;
            let n = ptr::read_unaligned(p);
            NodeView {
                hash: n.hash,
                coverage: n.coverage,
                degree: n.degree,
            }
        })
    }

    /// Iterate all edges (opaque iterator; yields EdgeView by value)
    pub fn iter_edges(&self) -> impl Iterator<Item = EdgeView> + '_ {
        let sz = size_of::<MappedEdge>();
        let total = self.edge_count;
        let mmap = &self.edge_mmap;
        (0..total).map(move |i| unsafe {
            let off = i * sz;
            let p = mmap.as_ptr().add(off) as *const MappedEdge;
            let e = ptr::read_unaligned(p);
            EdgeView {
                from_hash: e.from_hash,
                to_hash: e.to_hash,
                weight: e.weight,
                confidence: e.confidence,
            }
        })
    }

    /// Get statistics about memory usage
    pub fn get_memory_stats(&self) -> MemoryStats {
        let node_bytes = (self.node_count as u64) * (size_of::<MappedNode>() as u64);
        let edge_bytes = (self.edge_count as u64) * (size_of::<MappedEdge>() as u64);
        let index_bytes = self.index_mmap.len() as u64;

        let to_mb = |b: u64| b / (1024 * 1024);

        MemoryStats {
            nodes: self.node_count,
            edges: self.edge_count,
            node_memory_mb: to_mb(node_bytes),
            edge_memory_mb: to_mb(edge_bytes),
            index_memory_mb: to_mb(index_bytes),
            total_memory_mb: to_mb(node_bytes + edge_bytes + index_bytes),
        }
    }

    /// Flush all changes to disk
    pub fn flush(&mut self) -> Result<()> {
        self.node_mmap.flush()?;
        self.edge_mmap.flush()?;
        self.index_mmap.flush()?;
        Ok(())
    }

    /// Convert to standard AssemblyGraph (for compatibility)
    pub fn to_assembly_graph(&self) -> Result<AssemblyGraph> {
        let mut graph = AssemblyGraph::new();

        // Nodes
        for nv in self.iter_nodes() {
            let dummy_sequence = format!("NODE_{}", nv.hash);
            let kmer = CanonicalKmer::new(&dummy_sequence)
                .map_err(|e| anyhow::Error::msg(format!("Failed to create k-mer: {}", e)))?;
            let graph_node = GraphNode::new(kmer, 31); // k=31 (adjust as needed)
            graph.graph_fragment.add_node(graph_node);
        }

        // Edges
        for ev in self.iter_edges() {
            let graph_edge = GraphEdge::new(ev.from_hash, ev.to_hash, ev.weight as usize);
            graph.graph_fragment.add_edge(graph_edge);
        }

        Ok(graph)
    }

    // ------------------------- private helpers -------------------------

    fn create_sized_file(path: &str, size: usize) -> Result<File> {
        let file = OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(path)?;
        file.set_len(size as u64)?;
        Ok(file)
    }

    fn update_index(&mut self, hash: u64, offset: u64) -> Result<()> {
        let entry_sz = size_of::<IndexEntry>();
        let index_size = self.index_mmap.len() / entry_sz;
        if index_size == 0 {
            return Err(anyhow!("index file too small (0 slots)"));
        }
        let start = (hash as usize) % index_size;

        for i in 0..index_size {
            let pos = (start + i) % index_size;
            let entry_off = pos * entry_sz;
            unsafe {
                let p = self.index_mmap.as_mut_ptr().add(entry_off) as *mut IndexEntry;
                let mut entry = ptr::read_unaligned(p);
                if entry.hash == 0 || entry.hash == hash {
                    entry.hash = hash;
                    entry.offset = offset;
                    ptr::write_unaligned(p, entry);
                    return Ok(());
                }
            }
        }
        Err(anyhow!("Index full - cannot add more nodes"))
    }

    fn lookup_node_offset(&self, hash: u64) -> Option<u64> {
        let entry_sz = size_of::<IndexEntry>();
        let index_size = self.index_mmap.len() / entry_sz;
        if index_size == 0 {
            return None;
        }
        let start = (hash as usize) % index_size;

        for i in 0..index_size {
            let pos = (start + i) % index_size;
            let entry_off = pos * entry_sz;
            unsafe {
                let p = self.index_mmap.as_ptr().add(entry_off) as *const IndexEntry;
                let entry = ptr::read_unaligned(p);
                if entry.hash == hash {
                    return Some(entry.offset);
                }
                if entry.hash == 0 {
                    break; // empty slot => not found
                }
            }
        }
        None
    }
}

// --------------------------- Builder (public) ---------------------------

pub struct MappedGraphBuilder {
    base_path: String,
    estimated_nodes: usize,
    estimated_edges: usize,
}

impl MappedGraphBuilder {
    pub fn new<P: AsRef<Path>>(base_path: P) -> Self {
        Self {
            base_path: base_path.as_ref().to_string_lossy().to_string(),
            estimated_nodes: 1_000_000,
            estimated_edges: 2_000_000,
        }
    }

    pub fn with_capacity(mut self, nodes: usize, edges: usize) -> Self {
        self.estimated_nodes = nodes;
        self.estimated_edges = edges;
        self
    }

    pub fn build(self) -> Result<MappedGraph> {
        if let Some(parent) = Path::new(&self.base_path).parent() {
            std::fs::create_dir_all(parent)?;
        }
        MappedGraph::new(&self.base_path, self.estimated_nodes, self.estimated_edges)
    }

    /// Build from existing AssemblyGraph (conversion)
    pub fn build_from_graph(self, graph: &AssemblyGraph) -> Result<MappedGraph> {
        let node_count = graph.graph_fragment.nodes.len();
        let edge_count = graph.graph_fragment.edges.len();

        let mut mapped_graph = self
            .with_capacity(node_count * 2, edge_count * 2) // growth headroom
            .build()?;

        // Copy nodes (assuming graph.graph_fragment.nodes: HashMap<u64, GraphNodeLike>)
        for (hash, node) in &graph.graph_fragment.nodes {
            mapped_graph.add_node(*hash, node.coverage)?;
        }

        // Copy edges
        for edge in &graph.graph_fragment.edges {
            mapped_graph.add_edge(edge.from_hash, edge.to_hash, edge.weight)?;
        }

        mapped_graph.flush()?;
        Ok(mapped_graph)
    }
}

// -------------------------------- Tests --------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_memory_mapped_basic_operations() {
        let temp_dir = tempdir().unwrap();
        let mut graph = MappedGraph::new(temp_dir.path().join("test_graph"), 1000, 2000).unwrap();

        // Add some nodes
        graph.add_node(12345, 10).unwrap();
        graph.add_node(67890, 20).unwrap();

        // Add edges
        graph.add_edge(12345, 67890, 5).unwrap();

        // Test retrieval
        let node = graph.get_node(12345).unwrap();
        assert_eq!(node.hash, 12345);
        assert_eq!(node.coverage, 10);

        // Test iteration
        let nodes: Vec<NodeView> = graph.iter_nodes().collect();
        assert_eq!(nodes.len(), 2);

        let edges: Vec<EdgeView> = graph.iter_edges().collect();
        assert_eq!(edges.len(), 1);

        // Stats
        let stats = graph.get_memory_stats();
        assert_eq!(stats.nodes, 2);
        assert_eq!(stats.edges, 1);
        stats.print_summary();
    }

    #[test]
    fn test_large_graph_iteration() {
        let temp_dir = tempdir().unwrap();
        let mut graph =
            MappedGraph::new(temp_dir.path().join("large_test"), 10_000, 20_000).unwrap();

        for i in 0..1_000 {
            graph.add_node(i as u64, (i % 100) as u32).unwrap();
        }

        let count = graph.iter_nodes().count();
        assert_eq!(count, 1_000);
    }
}
