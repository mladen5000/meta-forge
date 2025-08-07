use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

/// Current screen in the TUI application
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Screen {
    MainMenu,
    FileSelection,
    Configuration,
    Analysis,
    Database,
    Results,
    Help,
    Error(String),
}

impl Default for Screen {
    fn default() -> Self {
        Screen::MainMenu
    }
}

/// Analysis configuration state
#[derive(Debug, Clone)]
pub struct AnalysisConfig {
    pub sample_name: String,
    pub input_files: Vec<PathBuf>,
    pub analysis_mode: String,
    pub k_mer_range: (usize, usize),
    pub min_coverage: u32,
    pub threads: Option<usize>,
    pub memory_limit: Option<usize>,
    pub output_dir: Option<PathBuf>,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            sample_name: String::new(),
            input_files: Vec::new(),
            analysis_mode: "Standard".to_string(),
            k_mer_range: (21, 127),
            min_coverage: 2,
            threads: None,
            memory_limit: None,
            output_dir: None,
        }
    }
}

/// Progress tracking for multiple operations
#[derive(Debug, Clone)]
pub struct OperationProgress {
    pub name: String,
    pub current: u64,
    pub total: Option<u64>,
    pub message: String,
    pub start_time: Instant,
    pub rate: Option<f64>,
    pub eta: Option<Duration>,
    pub status: OperationStatus,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OperationStatus {
    Pending,
    Running,
    Completed,
    Error(String),
    Cancelled,
}

impl OperationProgress {
    pub fn new(name: &str, total: Option<u64>) -> Self {
        Self {
            name: name.to_string(),
            current: 0,
            total,
            message: String::new(),
            start_time: Instant::now(),
            rate: None,
            eta: None,
            status: OperationStatus::Pending,
        }
    }

    pub fn update(&mut self, current: u64, message: Option<&str>) {
        self.current = current;
        if let Some(msg) = message {
            self.message = msg.to_string();
        }

        // Calculate rate and ETA
        let elapsed = self.start_time.elapsed().as_secs_f64();
        if elapsed > 0.0 {
            self.rate = Some(current as f64 / elapsed);

            if let (Some(total), Some(rate)) = (self.total, self.rate) {
                if rate > 0.0 && current < total {
                    let remaining = total - current;
                    let eta_seconds = remaining as f64 / rate;
                    self.eta = Some(Duration::from_secs_f64(eta_seconds));
                }
            }
        }
    }

    pub fn percentage(&self) -> Option<f64> {
        self.total.map(|total| {
            if total > 0 {
                (self.current as f64 / total as f64) * 100.0
            } else {
                0.0
            }
        })
    }
}

/// Database status information
#[derive(Debug, Clone)]
pub struct DatabaseInfo {
    pub path: Option<PathBuf>,
    pub initialized: bool,
    pub size: Option<u64>,
    pub tables: Vec<String>,
    pub last_updated: Option<String>,
}

impl Default for DatabaseInfo {
    fn default() -> Self {
        Self {
            path: None,
            initialized: false,
            size: None,
            tables: Vec::new(),
            last_updated: None,
        }
    }
}

/// File browser state
#[derive(Debug, Clone)]
pub struct FileBrowserState {
    pub current_dir: PathBuf,
    pub files: Vec<PathBuf>,
    pub selected_index: usize,
    pub selected_files: Vec<PathBuf>,
    pub filter: String,
    pub show_hidden: bool,
}

impl Default for FileBrowserState {
    fn default() -> Self {
        Self {
            current_dir: std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
            files: Vec::new(),
            selected_index: 0,
            selected_files: Vec::new(),
            filter: String::new(),
            show_hidden: false,
        }
    }
}

/// Analysis results data
#[derive(Debug, Clone, Default)]
pub struct ResultsData {
    pub sample_name: String,
    pub contigs_count: usize,
    pub total_length: u64,
    pub n50: u64,
    pub processing_time: Duration,
    pub assembly_stats: HashMap<String, String>,
    pub taxonomic_results: Vec<(String, f64)>,
    pub export_formats: Vec<String>,
}

/// Main application state
#[derive(Debug, Clone)]
pub struct AppState {
    /// Current screen being displayed
    pub current_screen: Screen,
    
    /// Previous screen for navigation
    pub previous_screen: Option<Screen>,
    
    /// Whether the application should quit
    pub should_quit: bool,
    
    /// Analysis configuration
    pub analysis_config: AnalysisConfig,
    
    /// Active operations and their progress
    pub operations: HashMap<String, OperationProgress>,
    
    /// Database status
    pub database_info: DatabaseInfo,
    
    /// File browser state
    pub file_browser: FileBrowserState,
    
    /// Analysis results
    pub results: ResultsData,
    
    /// Current error message if any
    pub error_message: Option<String>,
    
    /// Input state for forms
    pub input_mode: bool,
    
    /// Current input field name
    pub input_field: String,
    
    /// Temporary input buffer
    pub input_buffer: String,
    
    /// Status messages
    pub status_messages: Vec<String>,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            current_screen: Screen::default(),
            previous_screen: None,
            should_quit: false,
            analysis_config: AnalysisConfig::default(),
            operations: HashMap::new(),
            database_info: DatabaseInfo::default(),
            file_browser: FileBrowserState::default(),
            results: ResultsData::default(),
            error_message: None,
            input_mode: false,
            input_field: String::new(),
            input_buffer: String::new(),
            status_messages: Vec::new(),
        }
    }
}

impl AppState {
    pub fn new() -> Self {
        Self::default()
    }

    /// Navigate to a new screen
    pub fn navigate_to(&mut self, screen: Screen) {
        self.previous_screen = Some(self.current_screen.clone());
        self.current_screen = screen;
    }

    /// Go back to previous screen
    pub fn go_back(&mut self) {
        if let Some(previous) = self.previous_screen.take() {
            self.current_screen = previous;
        }
    }

    /// Add a status message
    pub fn add_status_message(&mut self, message: String) {
        self.status_messages.push(message);
        // Keep only last 100 messages
        if self.status_messages.len() > 100 {
            self.status_messages.remove(0);
        }
    }

    /// Set error state
    pub fn set_error(&mut self, error: String) {
        self.error_message = Some(error.clone());
        self.current_screen = Screen::Error(error);
    }

    /// Clear error state
    pub fn clear_error(&mut self) {
        self.error_message = None;
        if let Screen::Error(_) = self.current_screen {
            self.go_back();
        }
    }

    /// Start input mode for a field
    pub fn start_input(&mut self, field: &str) {
        self.input_mode = true;
        self.input_field = field.to_string();
        self.input_buffer.clear();
    }

    /// End input mode and return the input
    pub fn end_input(&mut self) -> String {
        self.input_mode = false;
        let input = self.input_buffer.clone();
        self.input_buffer.clear();
        self.input_field.clear();
        input
    }

    /// Add or update an operation
    pub fn add_operation(&mut self, name: &str, total: Option<u64>) {
        let progress = OperationProgress::new(name, total);
        self.operations.insert(name.to_string(), progress);
    }

    /// Update operation progress
    pub fn update_operation(&mut self, name: &str, current: u64, message: Option<&str>) {
        if let Some(operation) = self.operations.get_mut(name) {
            operation.update(current, message);
        }
    }

    /// Complete an operation
    pub fn complete_operation(&mut self, name: &str) {
        if let Some(operation) = self.operations.get_mut(name) {
            operation.status = OperationStatus::Completed;
        }
    }

    /// Set operation error
    pub fn error_operation(&mut self, name: &str, error: &str) {
        if let Some(operation) = self.operations.get_mut(name) {
            operation.status = OperationStatus::Error(error.to_string());
        }
    }

    /// Remove completed or errored operations
    pub fn cleanup_operations(&mut self) {
        self.operations.retain(|_, op| {
            !matches!(op.status, OperationStatus::Completed | OperationStatus::Error(_))
        });
    }
}