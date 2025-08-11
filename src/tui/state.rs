use std::collections::HashMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

/// Current screen in the TUI application
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum Screen {
    #[default]
    MainMenu,
    FileSelection,
    Configuration,
    Analysis,
    Database,
    Results,
    Help,
    Error(String),
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
#[derive(Debug, Clone, Default)]
pub struct DatabaseInfo {
    pub path: Option<PathBuf>,
    pub initialized: bool,
    pub size: Option<u64>,
    pub tables: Vec<String>,
    pub last_updated: Option<String>,
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
#[derive(Debug, Clone, Default)]
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
            !matches!(
                op.status,
                OperationStatus::Completed | OperationStatus::Error(_)
            )
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_screen_enum_transitions() {
        let mut state = AppState::new();

        // Test initial state
        assert_eq!(state.current_screen, Screen::MainMenu);
        assert!(state.previous_screen.is_none());

        // Test navigation
        state.navigate_to(Screen::FileSelection);
        assert_eq!(state.current_screen, Screen::FileSelection);
        assert_eq!(state.previous_screen, Some(Screen::MainMenu));

        // Test going back
        state.go_back();
        assert_eq!(state.current_screen, Screen::MainMenu);
    }

    #[test]
    fn test_analysis_config_default() {
        let config = AnalysisConfig::default();

        assert!(config.sample_name.is_empty());
        assert!(config.input_files.is_empty());
        assert_eq!(config.analysis_mode, "Standard");
        assert_eq!(config.k_mer_range, (21, 127));
        assert_eq!(config.min_coverage, 2);
        assert!(config.threads.is_none());
        assert!(config.memory_limit.is_none());
        assert!(config.output_dir.is_none());
    }

    #[test]
    fn test_operation_progress_new() {
        let progress = OperationProgress::new("test_operation", Some(100));

        assert_eq!(progress.name, "test_operation");
        assert_eq!(progress.current, 0);
        assert_eq!(progress.total, Some(100));
        assert!(progress.message.is_empty());
        assert_eq!(progress.status, OperationStatus::Pending);
        assert!(progress.rate.is_none());
        assert!(progress.eta.is_none());
    }

    #[test]
    fn test_operation_progress_update() {
        let mut progress = OperationProgress::new("test", Some(100));

        // Update progress
        progress.update(50, Some("Halfway done"));

        assert_eq!(progress.current, 50);
        assert_eq!(progress.message, "Halfway done");
        assert!(progress.rate.is_some());

        // Test percentage calculation
        let percentage = progress.percentage().unwrap();
        assert_eq!(percentage, 50.0);
    }

    #[test]
    fn test_operation_progress_without_total() {
        let mut progress = OperationProgress::new("indefinite", None);
        progress.update(42, Some("Processing..."));

        assert_eq!(progress.current, 42);
        assert!(progress.percentage().is_none());
    }

    #[test]
    fn test_database_info_default() {
        let db_info = DatabaseInfo::default();

        assert!(db_info.path.is_none());
        assert!(!db_info.initialized);
        assert!(db_info.size.is_none());
        assert!(db_info.tables.is_empty());
        assert!(db_info.last_updated.is_none());
    }

    #[test]
    fn test_file_browser_state_default() {
        let browser = FileBrowserState::default();

        assert!(browser.files.is_empty());
        assert_eq!(browser.selected_index, 0);
        assert!(browser.selected_files.is_empty());
        assert!(browser.filter.is_empty());
        assert!(!browser.show_hidden);
    }

    #[test]
    fn test_app_state_navigation() {
        let mut state = AppState::new();

        // Test navigation chain
        state.navigate_to(Screen::Configuration);
        state.navigate_to(Screen::Analysis);

        assert_eq!(state.current_screen, Screen::Analysis);
        assert_eq!(state.previous_screen, Some(Screen::Configuration));

        // Go back should return to Configuration, not MainMenu
        state.go_back();
        assert_eq!(state.current_screen, Screen::Configuration);
    }

    #[test]
    fn test_app_state_error_handling() {
        let mut state = AppState::new();

        state.set_error("Test error message".to_string());

        assert!(state.error_message.is_some());
        assert_eq!(state.error_message.as_ref().unwrap(), "Test error message");
        assert!(matches!(state.current_screen, Screen::Error(_)));

        // Clear error should restore previous screen
        state.previous_screen = Some(Screen::Configuration);
        state.clear_error();

        assert!(state.error_message.is_none());
        assert_eq!(state.current_screen, Screen::Configuration);
    }

    #[test]
    fn test_app_state_input_mode() {
        let mut state = AppState::new();

        assert!(!state.input_mode);
        assert!(state.input_field.is_empty());
        assert!(state.input_buffer.is_empty());

        // Start input
        state.start_input("sample_name");
        assert!(state.input_mode);
        assert_eq!(state.input_field, "sample_name");

        // Simulate typing
        state.input_buffer = "my_sample".to_string();

        // End input
        let input = state.end_input();
        assert!(!state.input_mode);
        assert!(state.input_field.is_empty());
        assert!(state.input_buffer.is_empty());
        assert_eq!(input, "my_sample");
    }

    #[test]
    fn test_status_messages() {
        let mut state = AppState::new();

        assert!(state.status_messages.is_empty());

        // Add some messages
        state.add_status_message("Started analysis".to_string());
        state.add_status_message("Processing reads".to_string());

        assert_eq!(state.status_messages.len(), 2);
        assert_eq!(state.status_messages[0], "Started analysis");
        assert_eq!(state.status_messages[1], "Processing reads");
    }

    #[test]
    fn test_status_messages_overflow() {
        let mut state = AppState::new();

        // Add more than 100 messages to test overflow behavior
        for i in 0..105 {
            state.add_status_message(format!("Message {}", i));
        }

        // Should keep only last 100 messages
        assert_eq!(state.status_messages.len(), 100);
        assert_eq!(state.status_messages[0], "Message 5"); // First 5 should be removed
        assert_eq!(state.status_messages[99], "Message 104");
    }

    #[test]
    fn test_operation_management() {
        let mut state = AppState::new();

        // Add operation
        state.add_operation("assembly", Some(1000));
        assert!(state.operations.contains_key("assembly"));

        // Update operation
        state.update_operation("assembly", 500, Some("Halfway"));
        let op = state.operations.get("assembly").unwrap();
        assert_eq!(op.current, 500);
        assert_eq!(op.message, "Halfway");

        // Complete operation
        state.complete_operation("assembly");
        let op = state.operations.get("assembly").unwrap();
        assert_eq!(op.status, OperationStatus::Completed);

        // Add error operation
        state.add_operation("failed_op", None);
        state.error_operation("failed_op", "Something went wrong");
        let op = state.operations.get("failed_op").unwrap();
        assert!(matches!(op.status, OperationStatus::Error(_)));

        // Cleanup should remove completed and errored operations
        let initial_count = state.operations.len();
        state.cleanup_operations();
        assert!(state.operations.len() < initial_count);
    }

    #[test]
    fn test_results_data_default() {
        let results = ResultsData::default();

        assert!(results.sample_name.is_empty());
        assert_eq!(results.contigs_count, 0);
        assert_eq!(results.total_length, 0);
        assert_eq!(results.n50, 0);
        assert_eq!(results.processing_time, Duration::ZERO);
        assert!(results.assembly_stats.is_empty());
        assert!(results.taxonomic_results.is_empty());
        assert!(results.export_formats.is_empty());
    }

    #[test]
    fn test_operation_status_equality() {
        assert_eq!(OperationStatus::Pending, OperationStatus::Pending);
        assert_eq!(OperationStatus::Running, OperationStatus::Running);
        assert_eq!(OperationStatus::Completed, OperationStatus::Completed);
        assert_eq!(
            OperationStatus::Error("test".to_string()),
            OperationStatus::Error("test".to_string())
        );
        assert_ne!(
            OperationStatus::Error("test1".to_string()),
            OperationStatus::Error("test2".to_string())
        );
    }
}
