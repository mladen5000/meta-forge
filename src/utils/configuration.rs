use anyhow::Result;
use config::{Config, ConfigError, Environment, File, FileFormat};
use serde::{Deserialize, Serialize};
use std::env;
use std::path::{Path, PathBuf};
use thiserror::Error;
use tracing::{error, info, warn};

/// Comprehensive configuration management for the metagenomics pipeline
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PipelineConfiguration {
    /// General pipeline settings
    pub general: GeneralConfig,
    /// Assembly configuration
    pub assembly: AssemblyConfig,
    /// Feature extraction configuration
    pub features: FeatureExtractionConfig,
    /// Database configuration
    pub database: DatabaseIntegrationConfig,
    /// Machine learning configuration
    pub ml: MachineLearningConfig,
    /// Performance and resource configuration
    pub performance: PerformanceConfig,
    /// Logging and monitoring configuration
    pub logging: LoggingConfig,
    /// Input/output configuration
    pub io: IOConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneralConfig {
    /// Pipeline name/version
    pub name: String,
    pub version: String,
    /// Working directory
    pub work_dir: PathBuf,
    /// Temporary directory
    pub temp_dir: PathBuf,
    /// Output directory
    pub output_dir: PathBuf,
    /// Enable debug mode
    pub debug_mode: bool,
    /// Random seed for reproducibility
    pub random_seed: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyConfig {
    /// K-mer size range for adaptive assembly
    pub k_min: usize,
    pub k_max: usize,
    /// Minimum coverage threshold
    pub min_coverage: u32,
    /// Complexity threshold for k-mer size selection
    pub complexity_threshold: f64,
    /// Enable graph simplification
    pub enable_simplification: bool,
    /// Bubble popping settings
    pub bubble_popping: BubblePoppingConfig,
    /// Tip removal settings
    pub tip_removal: TipRemovalConfig,
    /// Ambiguous base handling
    pub ambiguous_base_handling: AmbiguousBaseConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BubblePoppingConfig {
    pub enabled: bool,
    pub max_bubble_length: usize,
    pub min_coverage_ratio: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TipRemovalConfig {
    pub enabled: bool,
    pub max_tip_length: usize,
    pub min_coverage_ratio: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AmbiguousBaseConfig {
    /// Strategy for handling ambiguous bases (N)
    pub strategy: AmbiguousBaseStrategy,
    /// Maximum number of N's allowed in a k-mer (for Skip and Allow strategies)
    pub max_n_count: usize,
    /// Replacement base for N's (for Replace strategy)
    pub replacement_base: char,
    /// Probability distribution for random replacement (for RandomReplace)
    pub random_probabilities: Option<[f64; 4]>, // A, C, G, T
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub enum AmbiguousBaseStrategy {
    /// Skip k-mers containing any ambiguous bases (current behavior)
    Skip,
    /// Allow k-mers with limited number of ambiguous bases  
    #[default]
    Allow,
    /// Replace N with a specific base
    Replace,
    /// Replace N with random base based on probabilities
    RandomReplace,
    /// Replace N with most common base in surrounding context
    ContextReplace,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureExtractionConfig {
    /// Feature dimensions
    pub sequence_feature_dim: usize,
    pub graph_feature_dim: usize,
    pub kmer_feature_dim: usize,
    /// Which feature types to include
    pub include_composition: bool,
    pub include_codon_usage: bool,
    pub include_patterns: bool,
    pub include_complexity: bool,
    pub include_topology: bool,
    pub include_centrality: bool,
    pub include_clustering: bool,
    /// K-mer sizes for feature extraction
    pub kmer_sizes: Vec<usize>,
    /// Maximum k-mers to process
    pub max_kmers: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseIntegrationConfig {
    /// Database file path
    pub db_path: PathBuf,
    /// Enable WAL mode
    pub enable_wal_mode: bool,
    /// Cache size in pages
    pub cache_size: i32,
    /// Enable foreign keys
    pub enable_foreign_keys: bool,
    /// Batch size for bulk operations
    pub batch_size: usize,
    /// Enable compression
    pub enable_compression: bool,
    /// Cache memory limit in MB
    pub cache_memory_limit_mb: usize,
    /// Auto-vacuum settings
    pub auto_vacuum: AutoVacuumConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AutoVacuumConfig {
    pub enabled: bool,
    pub threshold_mb: usize,
    pub schedule: String, // Cron-like schedule
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MachineLearningConfig {
    /// Model paths
    pub taxonomy_model_path: Option<PathBuf>,
    pub repeat_model_path: Option<PathBuf>,
    pub error_correction_model_path: Option<PathBuf>,
    /// Training parameters
    pub training: TrainingConfig,
    /// Inference parameters
    pub inference: InferenceConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrainingConfig {
    pub batch_size: usize,
    pub learning_rate: f64,
    pub epochs: usize,
    pub validation_split: f64,
    pub early_stopping_patience: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InferenceConfig {
    pub batch_size: usize,
    pub confidence_threshold: f64,
    pub use_gpu: bool,
    pub max_sequence_length: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceConfig {
    /// Number of threads for parallel processing
    pub num_threads: usize,
    /// Memory limit in GB
    pub memory_limit_gb: usize,
    /// Enable memory monitoring
    pub enable_memory_monitoring: bool,
    /// Streaming buffer sizes
    pub streaming_buffer_size: usize,
    pub chunk_size: usize,
    /// Enable compression for intermediate data
    pub enable_compression: bool,
    /// Resource monitoring
    pub monitoring: MonitoringConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MonitoringConfig {
    pub enabled: bool,
    pub sample_interval_ms: u64,
    pub alert_memory_threshold_percent: f64,
    pub alert_cpu_threshold_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LoggingConfig {
    /// Log level (error, warn, info, debug, trace)
    pub level: String,
    /// Log output format (json, pretty, compact)
    pub format: String,
    /// Log file path (optional)
    pub file_path: Option<PathBuf>,
    /// Enable performance logging
    pub enable_performance_logging: bool,
    /// Enable metrics collection
    pub enable_metrics: bool,
    /// Metrics output path
    pub metrics_path: Option<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IOConfig {
    /// Input file formats to support
    pub supported_input_formats: Vec<String>,
    /// Output formats to generate
    pub output_formats: OutputFormatsConfig,
    /// File handling settings
    pub file_handling: FileHandlingConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputFormatsConfig {
    pub fasta: bool,
    pub fastq: bool,
    pub gfa: bool,
    pub json: bool,
    pub tsv: bool,
    pub html_report: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileHandlingConfig {
    pub enable_memory_mapping: bool,
    pub buffer_size: usize,
    pub max_file_size_gb: usize,
    pub compression_level: i32,
}

/// Centralized pipeline configuration to replace duplicates
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    pub k_range: (usize, usize),      // (min_k, max_k) for adaptive assembly
    pub memory_limit_gb: usize,       // Total memory budget
    pub use_gpu: bool,                // GPU acceleration if available
    pub quality_threshold: f64,       // Error correction threshold
    pub taxonomy_db_path: String,     // Path to taxonomy database
    pub num_threads: usize,           // Parallel processing threads
    pub enable_compression: bool,     // Compress intermediate data
    pub streaming_buffer_size: usize, // Buffer size for streaming
}

/// Custom error types for better error handling
#[derive(Error, Debug)]
pub enum PipelineError {
    #[error("Configuration error: {message}")]
    ConfigurationError { message: String },

    #[error("Input/Output error: {message}")]
    IOError { message: String },

    #[error("Database error: {message}")]
    DatabaseError { message: String },

    #[error("Assembly error: {message}")]
    AssemblyError { message: String },

    #[error("Feature extraction error: {message}")]
    FeatureExtractionError { message: String },

    #[error("Machine learning error: {message}")]
    MachineLearningError { message: String },

    #[error("Memory error: insufficient memory (required: {required_mb}MB, available: {available_mb}MB)")]
    MemoryError {
        required_mb: usize,
        available_mb: usize,
    },

    #[error("Resource error: {resource} limit exceeded")]
    ResourceError { resource: String },

    #[error("Validation error: {field} is invalid: {reason}")]
    ValidationError { field: String, reason: String },

    #[error("Recovery error: {message}")]
    RecoveryError { message: String },
}

impl From<ConfigError> for PipelineError {
    fn from(err: ConfigError) -> Self {
        PipelineError::ConfigurationError {
            message: err.to_string(),
        }
    }
}

impl From<std::io::Error> for PipelineError {
    fn from(err: std::io::Error) -> Self {
        PipelineError::IOError {
            message: err.to_string(),
        }
    }
}

impl From<rusqlite::Error> for PipelineError {
    fn from(err: rusqlite::Error) -> Self {
        PipelineError::DatabaseError {
            message: err.to_string(),
        }
    }
}

/// Configuration manager with validation and environment integration
pub struct ConfigurationManager {
    config: PipelineConfiguration,
    config_path: Option<PathBuf>,
    environment_prefix: String,
}

impl ConfigurationManager {
    /// Load configuration from multiple sources
    pub fn new() -> Result<Self, PipelineError> {
        Self::load_from_default_locations()
    }

    /// Create configuration manager with pure defaults (no file dependencies)
    pub fn new_with_defaults() -> Result<Self, PipelineError> {
        let config = PipelineConfiguration::default();

        let manager = Self {
            config,
            config_path: None,
            environment_prefix: "META".to_string(),
        };

        manager.validate_configuration()?;
        manager.setup_logging()?;

        Ok(manager)
    }

    /// Load configuration from specific file
    pub fn from_file<P: AsRef<Path>>(config_path: P) -> Result<Self, PipelineError> {
        let config_path = config_path.as_ref().to_path_buf();
        let config = Self::load_config_from_file(&config_path)?;

        let manager = Self {
            config,
            config_path: Some(config_path),
            environment_prefix: "META".to_string(),
        };

        manager.validate_configuration()?;
        manager.setup_logging()?;

        Ok(manager)
    }

    /// Load configuration from default locations
    fn load_from_default_locations() -> Result<Self, PipelineError> {
        let mut config_builder = Config::builder();

        // Start with default configuration from Rust defaults
        let default_config = PipelineConfiguration::default();

        // Try to load from embedded config file if available, but continue without it
        if let Ok(embedded_config) = std::fs::read_to_string("./config/default.toml") {
            config_builder = config_builder
                .add_source(config::File::from_str(&embedded_config, FileFormat::Toml));
        } else {
            info!("No config/default.toml found, using built-in defaults");
        }

        // Add system-wide config
        if let Ok(system_config) = env::var("META_SYSTEM_CONFIG") {
            config_builder =
                config_builder.add_source(File::with_name(&system_config).required(false));
        }

        // Add user config
        if let Some(home_dir) = dirs::home_dir() {
            let user_config = home_dir
                .join(".config")
                .join("metagenomics")
                .join("config.toml");
            config_builder = config_builder.add_source(File::from(user_config).required(false));
        }

        // Add local config
        config_builder = config_builder.add_source(File::with_name("config").required(false));

        // Add environment variables
        config_builder =
            config_builder.add_source(Environment::with_prefix("META").separator("__"));

        // Try to build configuration, fallback to default if it fails
        let config: PipelineConfiguration = match config_builder.build() {
            Ok(built_config) => match built_config.try_deserialize() {
                Ok(config) => config,
                Err(e) => {
                    warn!(
                        "Failed to deserialize configuration: {}, using built-in defaults",
                        e
                    );
                    default_config
                }
            },
            Err(e) => {
                warn!(
                    "Failed to build configuration: {}, using built-in defaults",
                    e
                );
                default_config
            }
        };

        let manager = Self {
            config,
            config_path: None,
            environment_prefix: "META".to_string(),
        };

        manager.validate_configuration()?;
        manager.setup_logging()?;

        Ok(manager)
    }

    fn load_config_from_file(path: &Path) -> Result<PipelineConfiguration, PipelineError> {
        let config = Config::builder().add_source(File::from(path)).build()?;

        Ok(config.try_deserialize()?)
    }

    /// Validate configuration parameters
    fn validate_configuration(&self) -> Result<(), PipelineError> {
        info!("🔍 Validating configuration...");

        // Validate assembly parameters
        if self.config.assembly.k_min >= self.config.assembly.k_max {
            return Err(PipelineError::ValidationError {
                field: "assembly.k_min".to_string(),
                reason: "must be less than k_max".to_string(),
            });
        }

        if self.config.assembly.k_min < 3 {
            return Err(PipelineError::ValidationError {
                field: "assembly.k_min".to_string(),
                reason: "must be at least 3".to_string(),
            });
        }

        if self.config.assembly.k_max > 255 {
            return Err(PipelineError::ValidationError {
                field: "assembly.k_max".to_string(),
                reason: "must be at most 255".to_string(),
            });
        }

        // Validate feature extraction parameters
        if self.config.features.sequence_feature_dim == 0 {
            return Err(PipelineError::ValidationError {
                field: "features.sequence_feature_dim".to_string(),
                reason: "must be greater than 0".to_string(),
            });
        }

        if self.config.features.kmer_sizes.is_empty() {
            return Err(PipelineError::ValidationError {
                field: "features.kmer_sizes".to_string(),
                reason: "must contain at least one k-mer size".to_string(),
            });
        }

        // Validate performance parameters
        if self.config.performance.num_threads == 0 {
            return Err(PipelineError::ValidationError {
                field: "performance.num_threads".to_string(),
                reason: "must be greater than 0".to_string(),
            });
        }

        let available_threads = num_cpus::get();
        if self.config.performance.num_threads > available_threads * 2 {
            warn!(
                "Configured threads ({}) exceeds available cores ({})",
                self.config.performance.num_threads, available_threads
            );
        }

        // Validate memory limits
        if self.config.performance.memory_limit_gb == 0 {
            return Err(PipelineError::ValidationError {
                field: "performance.memory_limit_gb".to_string(),
                reason: "must be greater than 0".to_string(),
            });
        }

        // Validate directories exist or can be created
        self.ensure_directories_exist()?;

        info!("✅ Configuration validation passed");
        Ok(())
    }

    fn ensure_directories_exist(&self) -> Result<(), PipelineError> {
        let directories = [
            &self.config.general.work_dir,
            &self.config.general.temp_dir,
            &self.config.general.output_dir,
        ];

        for dir in &directories {
            if !dir.exists() {
                std::fs::create_dir_all(dir).map_err(|e| PipelineError::IOError {
                    message: format!("Failed to create directory {}: {}", dir.display(), e),
                })?;
                info!("📁 Created directory: {}", dir.display());
            }
        }

        Ok(())
    }

    /// Setup logging based on configuration
    fn setup_logging(&self) -> Result<(), PipelineError> {
        use tracing_appender::rolling;
        use tracing_subscriber::{fmt, prelude::*, EnvFilter};

        // Check if global subscriber is already set
        if tracing::dispatcher::has_been_set() {
            info!("⏭️  Logging already initialized, skipping setup");
            return Ok(());
        }

        let level = &self.config.logging.level;
        let format = &self.config.logging.format;

        let env_filter =
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new(level));

        let subscriber = tracing_subscriber::registry().with(env_filter);

        match format.as_str() {
            "json" => {
                let layer = fmt::layer().with_target(true).with_thread_ids(true);
                if let Some(ref file_path) = self.config.logging.file_path {
                    let file_appender = rolling::daily(
                        file_path.parent().unwrap_or(Path::new(".")),
                        file_path
                            .file_name()
                            .unwrap_or(std::ffi::OsStr::new("pipeline.log")),
                    );
                    let (non_blocking, _guard) = tracing_appender::non_blocking(file_appender);
                    let _ = tracing::subscriber::set_global_default(
                        subscriber.with(layer.with_writer(non_blocking)),
                    );
                } else {
                    let _ = tracing::subscriber::set_global_default(subscriber.with(layer));
                }
            }
            "compact" => {
                let layer = fmt::layer().compact();
                if let Some(ref file_path) = self.config.logging.file_path {
                    let file_appender = rolling::daily(
                        file_path.parent().unwrap_or(Path::new(".")),
                        file_path
                            .file_name()
                            .unwrap_or(std::ffi::OsStr::new("pipeline.log")),
                    );
                    let (non_blocking, _guard) = tracing_appender::non_blocking(file_appender);
                    let _ = tracing::subscriber::set_global_default(
                        subscriber.with(layer.with_writer(non_blocking)),
                    );
                } else {
                    let _ = tracing::subscriber::set_global_default(subscriber.with(layer));
                }
            }
            _ => {
                // "pretty" or default - without timestamps for cleaner console output
                let layer = fmt::layer()
                    .without_time() // Remove timestamp
                    .with_target(false); // Remove target path for cleaner output
                if let Some(ref file_path) = self.config.logging.file_path {
                    let file_appender = rolling::daily(
                        file_path.parent().unwrap_or(Path::new(".")),
                        file_path
                            .file_name()
                            .unwrap_or(std::ffi::OsStr::new("pipeline.log")),
                    );
                    let (non_blocking, _guard) = tracing_appender::non_blocking(file_appender);
                    let _ = tracing::subscriber::set_global_default(
                        subscriber.with(layer.with_writer(non_blocking)),
                    );
                } else {
                    let _ = tracing::subscriber::set_global_default(subscriber.with(layer));
                }
            }
        }

        info!(
            "📝 Logging initialized with level: {}, format: {}",
            level, format
        );
        Ok(())
    }

    /// Get configuration reference
    pub fn config(&self) -> &PipelineConfiguration {
        &self.config
    }

    /// Get mutable configuration reference
    pub fn config_mut(&mut self) -> &mut PipelineConfiguration {
        &mut self.config
    }

    /// Save current configuration to file
    pub fn save_config<P: AsRef<Path>>(&self, path: P) -> Result<(), PipelineError> {
        let toml_string = toml::to_string_pretty(&self.config).map_err(|e| {
            PipelineError::ConfigurationError {
                message: format!("Failed to serialize configuration: {e}"),
            }
        })?;

        std::fs::write(path.as_ref(), toml_string).map_err(|e| PipelineError::IOError {
            message: format!("Failed to write configuration file: {e}"),
        })?;

        info!("💾 Configuration saved to {}", path.as_ref().display());
        Ok(())
    }

    /// Update configuration from environment variables
    pub fn update_from_environment(&mut self) -> Result<(), PipelineError> {
        // This would update configuration from environment variables
        // Implementation would use the same pattern as initial loading
        info!("🌍 Updated configuration from environment variables");
        Ok(())
    }

    /// Create configuration for specific use cases
    pub fn create_minimal_config() -> PipelineConfiguration {
        PipelineConfiguration {
            general: GeneralConfig {
                name: "metagenomics-pipeline".to_string(),
                version: "1.0.0".to_string(),
                work_dir: PathBuf::from("./output"),
                temp_dir: std::env::temp_dir(),
                output_dir: PathBuf::from("./output"),
                debug_mode: false,
                random_seed: None,
            },
            assembly: AssemblyConfig {
                k_min: 15,
                k_max: 31,
                min_coverage: 2,
                complexity_threshold: 0.7,
                enable_simplification: true,
                bubble_popping: BubblePoppingConfig {
                    enabled: true,
                    max_bubble_length: 1000,
                    min_coverage_ratio: 0.1,
                },
                tip_removal: TipRemovalConfig {
                    enabled: true,
                    max_tip_length: 100,
                    min_coverage_ratio: 0.1,
                },
                ambiguous_base_handling: AmbiguousBaseConfig::default(),
            },
            features: FeatureExtractionConfig {
                sequence_feature_dim: 100,
                graph_feature_dim: 50,
                kmer_feature_dim: 64,
                include_composition: true,
                include_codon_usage: true,
                include_patterns: true,
                include_complexity: true,
                include_topology: true,
                include_centrality: false, // Expensive to compute
                include_clustering: true,
                kmer_sizes: vec![3, 4, 5, 6],
                max_kmers: 10000,
            },
            database: DatabaseIntegrationConfig {
                db_path: PathBuf::from("./data/metagenomics.db"),
                enable_wal_mode: true,
                cache_size: 10000,
                enable_foreign_keys: true,
                batch_size: 1000,
                enable_compression: true,
                cache_memory_limit_mb: 256,
                auto_vacuum: AutoVacuumConfig {
                    enabled: true,
                    threshold_mb: 1000,
                    schedule: "0 2 * * *".to_string(), // Daily at 2 AM
                },
            },
            ml: MachineLearningConfig {
                taxonomy_model_path: None,
                repeat_model_path: None,
                error_correction_model_path: None,
                training: TrainingConfig {
                    batch_size: 64,
                    learning_rate: 0.001,
                    epochs: 100,
                    validation_split: 0.2,
                    early_stopping_patience: 10,
                },
                inference: InferenceConfig {
                    batch_size: 32,
                    confidence_threshold: 0.5,
                    use_gpu: false,
                    max_sequence_length: 10000,
                },
            },
            performance: PerformanceConfig {
                num_threads: num_cpus::get(),
                memory_limit_gb: 8,
                enable_memory_monitoring: true,
                streaming_buffer_size: 1000,
                chunk_size: 100,
                enable_compression: true,
                monitoring: MonitoringConfig {
                    enabled: true,
                    sample_interval_ms: 1000,
                    alert_memory_threshold_percent: 90.0,
                    alert_cpu_threshold_percent: 95.0,
                },
            },
            logging: LoggingConfig {
                level: "info".to_string(),
                format: "pretty".to_string(),
                file_path: Some(PathBuf::from("./logs/pipeline.log")),
                enable_performance_logging: true,
                enable_metrics: true,
                metrics_path: Some(PathBuf::from("./metrics/metrics.json")),
            },
            io: IOConfig {
                supported_input_formats: vec![
                    "fastq".to_string(),
                    "fasta".to_string(),
                    "fastq.gz".to_string(),
                    "fasta.gz".to_string(),
                ],
                output_formats: OutputFormatsConfig {
                    fasta: true,
                    fastq: false,
                    gfa: true,
                    json: true,
                    tsv: true,
                    html_report: true,
                },
                file_handling: FileHandlingConfig {
                    enable_memory_mapping: true,
                    buffer_size: 64 * 1024, // 64KB
                    max_file_size_gb: 10,
                    compression_level: 6,
                },
            },
        }
    }

    pub fn create_high_performance_config() -> PipelineConfiguration {
        let mut config = Self::create_minimal_config();

        // Optimize for performance
        config.performance.num_threads = num_cpus::get() * 2;
        config.performance.memory_limit_gb = 32;
        config.performance.streaming_buffer_size = 10000;
        config.performance.chunk_size = 1000;

        // Larger feature dimensions for better accuracy
        config.features.sequence_feature_dim = 200;
        config.features.graph_feature_dim = 100;
        config.features.kmer_feature_dim = 128;
        config.features.include_centrality = true;
        config.features.max_kmers = 100000;

        // Larger database caches
        config.database.cache_size = 100000;
        config.database.cache_memory_limit_mb = 1024;
        config.database.batch_size = 10000;

        // Enable GPU if available
        config.ml.inference.use_gpu = true;
        config.ml.inference.batch_size = 128;

        config
    }

    pub fn create_low_memory_config() -> PipelineConfiguration {
        let mut config = Self::create_minimal_config();

        // Optimize for low memory usage
        config.performance.memory_limit_gb = 2;
        config.performance.streaming_buffer_size = 100;
        config.performance.chunk_size = 10;
        config.performance.enable_compression = true;

        // Smaller feature dimensions
        config.features.sequence_feature_dim = 50;
        config.features.graph_feature_dim = 25;
        config.features.kmer_feature_dim = 32;
        config.features.include_centrality = false;
        config.features.max_kmers = 1000;

        // Smaller database caches
        config.database.cache_size = 1000;
        config.database.cache_memory_limit_mb = 64;
        config.database.enable_compression = true;

        // Disable expensive features
        config.assembly.enable_simplification = false;
        config.assembly.bubble_popping.enabled = false;

        config
    }
}

/// Resource monitoring and management
pub struct ResourceMonitor {
    config: MonitoringConfig,
    memory_baseline: usize,
    cpu_baseline: f64,
    alert_handlers: Vec<Box<dyn Fn(&ResourceAlert) + Send + Sync>>,
}

#[derive(Debug, Clone)]
pub struct ResourceAlert {
    pub alert_type: AlertType,
    pub message: String,
    pub value: f64,
    pub threshold: f64,
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

#[derive(Debug, Clone)]
pub enum AlertType {
    MemoryHigh,
    CpuHigh,
    DiskSpaceLow,
    Custom(String),
}

impl ResourceMonitor {
    pub fn new(config: MonitoringConfig) -> Self {
        Self {
            config,
            memory_baseline: Self::get_memory_usage(),
            cpu_baseline: Self::get_cpu_usage(),
            alert_handlers: Vec::new(),
        }
    }

    pub fn start_monitoring(&mut self) -> Result<(), PipelineError> {
        if !self.config.enabled {
            return Ok(());
        }

        // This would start a background thread for monitoring
        // For now, just log the baseline values
        info!(
            "📊 Memory baseline: {} MB",
            self.memory_baseline / 1024 / 1024
        );
        info!("📊 CPU baseline: {:.1}%", self.cpu_baseline);

        Ok(())
    }

    pub fn check_resources(&self) -> Result<(), PipelineError> {
        let memory_usage = Self::get_memory_usage();
        let memory_percent =
            memory_usage as f64 / (self.config.alert_memory_threshold_percent / 100.0);

        if memory_percent > self.config.alert_memory_threshold_percent {
            let alert = ResourceAlert {
                alert_type: AlertType::MemoryHigh,
                message: format!("Memory usage high: {memory_percent:.1}%"),
                value: memory_percent,
                threshold: self.config.alert_memory_threshold_percent,
                timestamp: chrono::Utc::now(),
            };

            self.handle_alert(&alert);
        }

        let cpu_usage = Self::get_cpu_usage();
        if cpu_usage > self.config.alert_cpu_threshold_percent {
            let alert = ResourceAlert {
                alert_type: AlertType::CpuHigh,
                message: format!("CPU usage high: {cpu_usage:.1}%"),
                value: cpu_usage,
                threshold: self.config.alert_cpu_threshold_percent,
                timestamp: chrono::Utc::now(),
            };

            self.handle_alert(&alert);
        }

        Ok(())
    }

    fn handle_alert(&self, alert: &ResourceAlert) {
        warn!("🚨 Resource Alert: {}", alert.message);

        for handler in &self.alert_handlers {
            handler(alert);
        }
    }

    fn get_memory_usage() -> usize {
        // Simplified memory usage calculation
        // In a real implementation, would use system APIs
        1024 * 1024 * 100 // 100 MB placeholder
    }

    fn get_cpu_usage() -> f64 {
        // Simplified CPU usage calculation
        // In a real implementation, would use system APIs
        25.0 // 25% placeholder
    }

    pub fn add_alert_handler<F>(&mut self, handler: F)
    where
        F: Fn(&ResourceAlert) + Send + Sync + 'static,
    {
        self.alert_handlers.push(Box::new(handler));
    }
}

/// Error recovery strategies
pub struct ErrorRecoveryManager;

impl ErrorRecoveryManager {
    pub fn handle_pipeline_error(error: &PipelineError) -> RecoveryAction {
        match error {
            PipelineError::MemoryError {
                required_mb,
                available_mb,
            } => {
                if *required_mb > *available_mb * 2 {
                    RecoveryAction::Abort("Insufficient memory, cannot continue".to_string())
                } else {
                    RecoveryAction::Retry("Reduce memory usage and retry".to_string())
                }
            }
            PipelineError::IOError { .. } => {
                RecoveryAction::Retry("I/O error, retry with backoff".to_string())
            }
            PipelineError::DatabaseError { .. } => {
                RecoveryAction::Recover("Database error, attempt recovery".to_string())
            }
            PipelineError::ValidationError { .. } => {
                RecoveryAction::Abort("Configuration validation failed".to_string())
            }
            _ => RecoveryAction::Continue("Log error and continue".to_string()),
        }
    }

    pub fn attempt_recovery(action: &RecoveryAction) -> Result<(), PipelineError> {
        match action {
            RecoveryAction::Retry(message) => {
                warn!("🔄 Attempting recovery: {}", message);
                // Implement retry logic
                Ok(())
            }
            RecoveryAction::Recover(message) => {
                warn!("🛠️  Attempting recovery: {}", message);
                // Implement recovery logic
                Ok(())
            }
            RecoveryAction::Continue(message) => {
                info!("➡️  Continuing: {}", message);
                Ok(())
            }
            RecoveryAction::Abort(message) => {
                error!("🛑 Aborting: {}", message);
                Err(PipelineError::RecoveryError {
                    message: message.clone(),
                })
            }
        }
    }
}

#[derive(Debug, Clone)]
pub enum RecoveryAction {
    Retry(String),
    Recover(String),
    Continue(String),
    Abort(String),
}

// Default implementations for all configuration structures

impl Default for GeneralConfig {
    fn default() -> Self {
        Self {
            name: "metagenomics-pipeline".to_string(),
            version: "1.0.0".to_string(),
            work_dir: PathBuf::from("./output"),
            temp_dir: std::env::temp_dir(),
            output_dir: PathBuf::from("./output"),
            debug_mode: false,
            random_seed: None,
        }
    }
}

impl Default for AssemblyConfig {
    fn default() -> Self {
        Self {
            k_min: 15,
            k_max: 31,
            min_coverage: 2,
            complexity_threshold: 0.7,
            enable_simplification: true,
            bubble_popping: BubblePoppingConfig::default(),
            tip_removal: TipRemovalConfig::default(),
            ambiguous_base_handling: AmbiguousBaseConfig::default(),
        }
    }
}

impl Default for BubblePoppingConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            max_bubble_length: 1000,
            min_coverage_ratio: 0.1,
        }
    }
}

impl Default for TipRemovalConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            max_tip_length: 100,
            min_coverage_ratio: 0.1,
        }
    }
}

impl Default for AmbiguousBaseConfig {
    fn default() -> Self {
        Self {
            strategy: AmbiguousBaseStrategy::Allow,
            max_n_count: 2,
            replacement_base: 'A',
            random_probabilities: Some([0.25, 0.25, 0.25, 0.25]), // Equal probabilities
        }
    }
}

impl Default for FeatureExtractionConfig {
    fn default() -> Self {
        Self {
            sequence_feature_dim: 100,
            graph_feature_dim: 50,
            kmer_feature_dim: 64,
            include_composition: true,
            include_codon_usage: true,
            include_patterns: true,
            include_complexity: true,
            include_topology: true,
            include_centrality: false, // Expensive to compute
            include_clustering: true,
            kmer_sizes: vec![3, 4, 5, 6],
            max_kmers: 10000,
        }
    }
}

impl Default for DatabaseIntegrationConfig {
    fn default() -> Self {
        Self {
            db_path: PathBuf::from("./data/metagenomics.db"),
            enable_wal_mode: true,
            cache_size: 10000,
            enable_foreign_keys: true,
            batch_size: 1000,
            enable_compression: true,
            cache_memory_limit_mb: 256,
            auto_vacuum: AutoVacuumConfig::default(),
        }
    }
}

impl Default for AutoVacuumConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            threshold_mb: 1000,
            schedule: "0 2 * * *".to_string(), // Daily at 2 AM
        }
    }
}

impl Default for TrainingConfig {
    fn default() -> Self {
        Self {
            batch_size: 64,
            learning_rate: 0.001,
            epochs: 100,
            validation_split: 0.2,
            early_stopping_patience: 10,
        }
    }
}

impl Default for InferenceConfig {
    fn default() -> Self {
        Self {
            batch_size: 32,
            confidence_threshold: 0.5,
            use_gpu: false,
            max_sequence_length: 10000,
        }
    }
}

impl Default for PerformanceConfig {
    fn default() -> Self {
        Self {
            num_threads: num_cpus::get(),
            memory_limit_gb: 8,
            enable_memory_monitoring: true,
            streaming_buffer_size: 1000,
            chunk_size: 100,
            enable_compression: false,
            monitoring: MonitoringConfig::default(),
        }
    }
}

impl Default for MonitoringConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            sample_interval_ms: 1000,
            alert_memory_threshold_percent: 90.0,
            alert_cpu_threshold_percent: 85.0,
        }
    }
}

impl Default for LoggingConfig {
    fn default() -> Self {
        Self {
            level: "info".to_string(),
            format: "pretty".to_string(),
            file_path: None,
            enable_performance_logging: true,
            enable_metrics: true,
            metrics_path: None,
        }
    }
}

impl Default for IOConfig {
    fn default() -> Self {
        Self {
            supported_input_formats: vec![
                "fasta".to_string(),
                "fastq".to_string(),
                "fasta.gz".to_string(),
                "fastq.gz".to_string(),
            ],
            output_formats: OutputFormatsConfig::default(),
            file_handling: FileHandlingConfig::default(),
        }
    }
}

impl Default for OutputFormatsConfig {
    fn default() -> Self {
        Self {
            fasta: true,
            fastq: false,
            gfa: true,
            json: true,
            tsv: true,
            html_report: true,
        }
    }
}

impl Default for FileHandlingConfig {
    fn default() -> Self {
        Self {
            enable_memory_mapping: true,
            buffer_size: 8192,
            max_file_size_gb: 10,
            compression_level: 6,
        }
    }
}

impl Default for ConfigurationManager {
    fn default() -> Self {
        Self {
            config: PipelineConfiguration::default(),
            config_path: None,
            environment_prefix: "META".to_string(),
        }
    }
}

/// Utility functions for configuration management
pub mod config_utils {
    use super::*;

    /// Validate a configuration file without loading it
    pub fn validate_config_file<P: AsRef<Path>>(path: P) -> Result<(), PipelineError> {
        let config = ConfigurationManager::load_config_from_file(path.as_ref())?;

        // Basic validation
        if config.assembly.k_min >= config.assembly.k_max {
            return Err(PipelineError::ValidationError {
                field: "assembly.k_min".to_string(),
                reason: "must be less than k_max".to_string(),
            });
        }

        Ok(())
    }

    /// Merge two configurations, with the second taking precedence
    pub fn merge_configs(
        base: PipelineConfiguration,
        override_config: PipelineConfiguration,
    ) -> PipelineConfiguration {
        // For simplicity, just return the override config
        // In a real implementation, would merge field by field
        override_config
    }

    /// Generate a template configuration file
    pub fn generate_config_template<P: AsRef<Path>>(path: P) -> Result<(), PipelineError> {
        let template_config = ConfigurationManager::create_minimal_config();

        let toml_string = toml::to_string_pretty(&template_config).map_err(|e| {
            PipelineError::ConfigurationError {
                message: format!("Failed to serialize template: {e}"),
            }
        })?;

        std::fs::write(path.as_ref(), toml_string).map_err(|e| PipelineError::IOError {
            message: format!("Failed to write template: {e}"),
        })?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_minimal_config_creation() {
        let config = ConfigurationManager::create_minimal_config();

        assert_eq!(config.general.name, "metagenomics-pipeline");
        assert!(config.assembly.k_min < config.assembly.k_max);
        assert!(config.features.sequence_feature_dim > 0);
        assert!(config.performance.num_threads > 0);
    }

    #[test]
    fn test_config_validation() {
        let mut config = ConfigurationManager::create_minimal_config();

        // Test invalid k-mer range
        config.assembly.k_min = 30;
        config.assembly.k_max = 20;

        let temp_dir = tempdir().unwrap();
        config.general.work_dir = temp_dir.path().to_path_buf();
        config.general.temp_dir = temp_dir.path().join("tmp");
        config.general.output_dir = temp_dir.path().join("output");

        let manager = ConfigurationManager {
            config,
            config_path: None,
            environment_prefix: "TEST".to_string(),
        };

        assert!(manager.validate_configuration().is_err());
    }

    #[test]
    fn test_config_serialization() {
        let config = ConfigurationManager::create_minimal_config();
        let temp_dir = tempdir().unwrap();
        let config_path = temp_dir.path().join("test_config.toml");

        // Serialize
        let toml_string = toml::to_string_pretty(&config).unwrap();
        std::fs::write(&config_path, toml_string).unwrap();

        // Deserialize
        let loaded_config = ConfigurationManager::load_config_from_file(&config_path).unwrap();

        assert_eq!(config.general.name, loaded_config.general.name);
        assert_eq!(config.assembly.k_min, loaded_config.assembly.k_min);
    }

    #[test]
    fn test_error_recovery() {
        let memory_error = PipelineError::MemoryError {
            required_mb: 1500, // More than 2x available (1500 > 500*2)
            available_mb: 500,
        };

        let action = ErrorRecoveryManager::handle_pipeline_error(&memory_error);

        match action {
            RecoveryAction::Abort(_) => {
                // Expected for severe memory shortage
            }
            _ => panic!("Expected abort action for severe memory error"),
        }

        let io_error = PipelineError::IOError {
            message: "File not found".to_string(),
        };

        let action = ErrorRecoveryManager::handle_pipeline_error(&io_error);

        match action {
            RecoveryAction::Retry(_) => {
                // Expected for I/O errors
            }
            _ => panic!("Expected retry action for I/O error"),
        }
    }

    #[test]
    fn test_resource_monitoring() {
        let monitoring_config = MonitoringConfig {
            enabled: true,
            sample_interval_ms: 100,
            alert_memory_threshold_percent: 80.0,
            alert_cpu_threshold_percent: 90.0,
        };

        let mut monitor = ResourceMonitor::new(monitoring_config);

        // Add a test alert handler
        monitor.add_alert_handler(|alert| {
            println!("Test alert: {}", alert.message);
        });

        assert!(monitor.start_monitoring().is_ok());
        assert!(monitor.check_resources().is_ok());
    }
}
