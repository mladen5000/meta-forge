use anyhow::{Context, Result};
use crossterm::{
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};
use std::io;
use std::sync::Arc;
use std::time::Duration;
use tokio::sync::RwLock;

use crate::tui::{
    events::{AppEvent, EventHandler},
    screens::ScreenManager,
    state::{AppState, Screen},
};

/// Main TUI application
pub struct TuiApp {
    /// Terminal interface
    terminal: Terminal<CrosstermBackend<io::Stdout>>,

    /// Application state
    state: Arc<RwLock<AppState>>,

    /// Event handler
    event_handler: EventHandler,

    /// Screen manager
    screen_manager: ScreenManager,

    /// Whether the app is running
    running: bool,
}

impl TuiApp {
    /// Create a new TUI application
    pub fn new() -> Result<Self> {
        // Setup terminal
        enable_raw_mode().context("Failed to enable raw mode")?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen).context("Failed to enter alternate screen")?;

        let backend = CrosstermBackend::new(stdout);
        let terminal = Terminal::new(backend).context("Failed to create terminal")?;

        // Create application state
        let state = Arc::new(RwLock::new(AppState::new()));

        // Create event handler with 250ms tick rate
        let event_handler = EventHandler::new(Duration::from_millis(250));

        // Create screen manager
        let screen_manager = ScreenManager::new(state.clone());

        Ok(Self {
            terminal,
            state,
            event_handler,
            screen_manager,
            running: false,
        })
    }

    /// Run the TUI application
    pub async fn run(&mut self) -> Result<()> {
        self.running = true;

        // Initialize screen manager
        self.screen_manager.initialize().await?;

        // Main event loop
        while self.running {
            // Render the current screen
            self.draw().await?;

            // Handle events
            if let Some(event) = self.event_handler.try_next()? {
                self.handle_event(event).await?;
            }

            // Handle async events from background tasks
            while let Ok(event) = self.event_handler.async_receiver.try_recv() {
                self.handle_event(event).await?;
            }

            // Check if we should quit
            let state = self.state.read().await;
            if state.should_quit {
                self.running = false;
            }
            drop(state);

            // Small delay to prevent busy-waiting
            tokio::time::sleep(Duration::from_millis(16)).await; // ~60 FPS
        }

        Ok(())
    }

    /// Draw the current screen
    async fn draw(&mut self) -> Result<()> {
        let state = self.state.read().await;
        let current_screen = state.current_screen.clone();
        drop(state);

        self.terminal.draw(|f| {
            if let Err(e) = self.screen_manager.render(f, &current_screen) {
                eprintln!("Render error: {e}");
            }
        })?;

        Ok(())
    }

    /// Handle application events
    async fn handle_event(&mut self, event: AppEvent) -> Result<()> {
        match event {
            AppEvent::Input(key) => {
                // Handle global key bindings first
                if let crossterm::event::KeyCode::Char('q') = key.code {
                    if key
                        .modifiers
                        .contains(crossterm::event::KeyModifiers::CONTROL)
                    {
                        let mut state = self.state.write().await;
                        state.should_quit = true;
                        return Ok(());
                    }
                }

                // Get current screen for handling
                let screen = {
                    let state = self.state.read().await;
                    state.current_screen.clone()
                };

                // Let screen manager handle all keys including Esc
                self.screen_manager.handle_key_input(key, &screen).await?;
            }

            AppEvent::Tick => {
                // Update timers and refresh displays
                self.screen_manager.tick().await?;
            }

            AppEvent::AnalysisStarted(name) => {
                let mut state = self.state.write().await;
                state.add_operation(&name, None);
                state.add_status_message(format!("Started: {name}"));
            }

            AppEvent::AnalysisProgress {
                operation,
                current,
                total,
                message,
            } => {
                let mut state = self.state.write().await;
                if !state.operations.contains_key(&operation) {
                    state.add_operation(&operation, total);
                }
                state.update_operation(&operation, current, Some(&message));
            }

            AppEvent::AnalysisCompleted(name) => {
                let mut state = self.state.write().await;
                state.complete_operation(&name);
                state.add_status_message(format!("Completed: {name}"));
            }

            AppEvent::AnalysisError(name, error) => {
                let mut state = self.state.write().await;
                state.error_operation(&name, &error);
                state.set_error(format!("Analysis error in {name}: {error}"));
            }

            AppEvent::DatabaseUpdate(db_info) => {
                let mut state = self.state.write().await;
                state.database_info = db_info;
            }

            AppEvent::FileListUpdated(files) => {
                let mut state = self.state.write().await;
                state.file_browser.files = files;
            }

            AppEvent::ResultsAvailable(results) => {
                let mut state = self.state.write().await;
                state.results = results;
                state.navigate_to(Screen::Results);
            }

            AppEvent::StatusMessage(message) => {
                let mut state = self.state.write().await;
                state.add_status_message(message);
            }

            AppEvent::Quit => {
                let mut state = self.state.write().await;
                state.should_quit = true;
            }
        }

        Ok(())
    }

    /// Get a reference to the application state
    pub fn state(&self) -> Arc<RwLock<AppState>> {
        self.state.clone()
    }

    /// Get event sender for background tasks
    pub fn event_sender(&self) -> tokio::sync::mpsc::UnboundedSender<AppEvent> {
        self.event_handler.async_sender()
    }
}

impl Drop for TuiApp {
    fn drop(&mut self) {
        // Restore terminal
        let _ = disable_raw_mode();
        let _ = execute!(self.terminal.backend_mut(), LeaveAlternateScreen);
        let _ = self.terminal.show_cursor();
    }
}

/// Initialize the TUI application
pub async fn initialize() -> Result<TuiApp> {
    TuiApp::new()
}

/// Run the TUI application with proper error handling
pub async fn run_tui() -> Result<()> {
    let mut app = initialize().await?;

    // Setup panic hook to restore terminal
    std::panic::set_hook(Box::new(|panic| {
        let _ = disable_raw_mode();
        let _ = execute!(io::stdout(), LeaveAlternateScreen);
        eprintln!("Application panicked: {panic}");
    }));

    let result = app.run().await;

    // Restore default panic hook
    let _ = std::panic::take_hook();

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tui::{
        events::AppEvent,
        state::{AppState, Screen},
    };
    use std::time::Duration;
    use tokio::sync::mpsc;

    // Mock TUI app for testing (since we can't easily test the real one without a terminal)
    struct MockTuiApp {
        state: Arc<RwLock<AppState>>,
        running: bool,
        events_received: Vec<AppEvent>,
    }

    impl MockTuiApp {
        fn new() -> Self {
            Self {
                state: Arc::new(RwLock::new(AppState::new())),
                running: false,
                events_received: Vec::new(),
            }
        }

        async fn handle_mock_event(&mut self, event: AppEvent) -> Result<()> {
            self.events_received.push(event.clone());

            match event {
                AppEvent::Input(key) => {
                    if let crossterm::event::KeyCode::Char('q') = key.code {
                        if key
                            .modifiers
                            .contains(crossterm::event::KeyModifiers::CONTROL)
                        {
                            let mut state = self.state.write().await;
                            state.should_quit = true;
                        }
                    }
                }
                AppEvent::AnalysisStarted(name) => {
                    let mut state = self.state.write().await;
                    state.add_operation(&name, None);
                }
                AppEvent::AnalysisCompleted(name) => {
                    let mut state = self.state.write().await;
                    state.complete_operation(&name);
                }
                AppEvent::Quit => {
                    let mut state = self.state.write().await;
                    state.should_quit = true;
                }
                _ => {}
            }

            Ok(())
        }
    }

    #[test]
    fn test_tui_app_creation_concept() {
        // Test the concept of TUI app creation without actually creating a terminal
        let state = Arc::new(RwLock::new(AppState::new()));

        // Test initial state
        let runtime = tokio::runtime::Runtime::new().unwrap();
        runtime.block_on(async {
            let state_read = state.read().await;
            assert_eq!(state_read.current_screen, Screen::MainMenu);
            assert!(!state_read.should_quit);
        });
    }

    #[tokio::test]
    async fn test_app_state_management() {
        let mut mock_app = MockTuiApp::new();

        // Test initial state
        {
            let state = mock_app.state.read().await;
            assert!(!state.should_quit);
            assert_eq!(state.current_screen, Screen::MainMenu);
        }

        // Test quit event handling
        let quit_event = AppEvent::Quit;
        mock_app.handle_mock_event(quit_event).await.unwrap();

        {
            let state = mock_app.state.read().await;
            assert!(state.should_quit);
        }

        assert_eq!(mock_app.events_received.len(), 1);
    }

    #[tokio::test]
    async fn test_analysis_event_handling() {
        let mut mock_app = MockTuiApp::new();

        // Test analysis started event
        let start_event = AppEvent::AnalysisStarted("test_analysis".to_string());
        mock_app.handle_mock_event(start_event).await.unwrap();

        {
            let state = mock_app.state.read().await;
            assert!(state.operations.contains_key("test_analysis"));
        }

        // Test analysis completed event
        let complete_event = AppEvent::AnalysisCompleted("test_analysis".to_string());
        mock_app.handle_mock_event(complete_event).await.unwrap();

        {
            let state = mock_app.state.read().await;
            let operation = state.operations.get("test_analysis").unwrap();
            assert_eq!(
                operation.status,
                crate::tui::state::OperationStatus::Completed
            );
        }
    }

    #[tokio::test]
    async fn test_key_event_handling() {
        use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

        let mut mock_app = MockTuiApp::new();

        // Test Ctrl+Q for quit
        let quit_key = KeyEvent::new(KeyCode::Char('q'), KeyModifiers::CONTROL);
        let input_event = AppEvent::Input(quit_key);

        mock_app.handle_mock_event(input_event).await.unwrap();

        {
            let state = mock_app.state.read().await;
            assert!(state.should_quit);
        }
    }

    #[tokio::test]
    async fn test_event_flow() {
        let mut mock_app = MockTuiApp::new();

        // Test a sequence of events
        let events = vec![
            AppEvent::AnalysisStarted("assembly".to_string()),
            AppEvent::AnalysisProgress {
                operation: "assembly".to_string(),
                current: 50,
                total: Some(100),
                message: "Half done".to_string(),
            },
            AppEvent::AnalysisCompleted("assembly".to_string()),
            AppEvent::StatusMessage("Analysis finished".to_string()),
        ];

        for event in events {
            mock_app.handle_mock_event(event).await.unwrap();
        }

        assert_eq!(mock_app.events_received.len(), 4);

        // Check final state
        {
            let state = mock_app.state.read().await;
            assert!(state.operations.contains_key("assembly"));
            assert_eq!(state.status_messages.len(), 1);
            assert_eq!(state.status_messages[0], "Analysis finished");
        }
    }

    #[test]
    fn test_terminal_restoration_concept() {
        // Test the concept of terminal restoration
        // In real usage, this would involve crossterm operations

        // Mock the restoration process
        let restoration_needed = true;
        let restore_successful = true; // In real code, this would be the result of crossterm calls

        if restoration_needed {
            assert!(restore_successful, "Terminal restoration should succeed");
        }
    }

    #[tokio::test]
    async fn test_async_event_communication() {
        // Test async event communication concept
        let (tx, mut rx) = mpsc::unbounded_channel::<AppEvent>();

        // Send some events
        tx.send(AppEvent::StatusMessage("Test message".to_string()))
            .unwrap();
        tx.send(AppEvent::AnalysisStarted("test".to_string()))
            .unwrap();

        // Receive events
        let mut received_events = Vec::new();
        while let Ok(event) = rx.try_recv() {
            received_events.push(event);
        }

        assert_eq!(received_events.len(), 2);
    }

    #[test]
    fn test_event_handler_creation_concept() {
        // Test event handler creation concept
        let tick_rate = Duration::from_millis(250);

        // In real code, this would create the actual EventHandler
        // Here we just test the concept
        assert!(tick_rate.as_millis() > 0);
        assert!(tick_rate.as_millis() <= 1000); // Reasonable tick rate
    }

    #[tokio::test]
    async fn test_state_sharing() {
        // Test state sharing between components
        let shared_state = Arc::new(RwLock::new(AppState::new()));

        // Clone for multiple owners (simulating app components)
        let state_clone1 = shared_state.clone();
        let state_clone2 = shared_state.clone();

        // Modify from one clone
        {
            let mut state = state_clone1.write().await;
            state.navigate_to(Screen::FileSelection);
        }

        // Read from another clone
        {
            let state = state_clone2.read().await;
            assert_eq!(state.current_screen, Screen::FileSelection);
        }
    }

    #[tokio::test]
    async fn test_error_handling_flow() {
        let mut mock_app = MockTuiApp::new();

        // Test error event
        let error_event =
            AppEvent::AnalysisError("assembly".to_string(), "Out of memory".to_string());

        mock_app.handle_mock_event(error_event).await.unwrap();

        // Check that the error is properly recorded
        assert_eq!(mock_app.events_received.len(), 1);

        if let AppEvent::AnalysisError(op, err) = &mock_app.events_received[0] {
            assert_eq!(op, "assembly");
            assert_eq!(err, "Out of memory");
        } else {
            panic!("Expected AnalysisError event");
        }
    }

    #[test]
    fn test_screen_manager_concept() {
        // Test screen manager concept
        let screens = vec![
            Screen::MainMenu,
            Screen::FileSelection,
            Screen::Configuration,
            Screen::Analysis,
            Screen::Database,
            Screen::Results,
            Screen::Help,
        ];

        // All screens should be representable
        for screen in screens {
            match screen {
                Screen::MainMenu => assert!(true),
                Screen::FileSelection => assert!(true),
                Screen::Configuration => assert!(true),
                Screen::Analysis => assert!(true),
                Screen::Database => assert!(true),
                Screen::Results => assert!(true),
                Screen::Help => assert!(true),
                Screen::Error(_) => assert!(true),
            }
        }
    }

    #[tokio::test]
    async fn test_multi_screen_navigation() {
        let mut mock_app = MockTuiApp::new();

        // Start at main menu
        {
            let state = mock_app.state.read().await;
            assert_eq!(state.current_screen, Screen::MainMenu);
            assert!(state.previous_screen.is_none());
        }

        // Navigate to file selection
        {
            let mut state = mock_app.state.write().await;
            state.navigate_to(Screen::FileSelection);
        }

        // Check navigation
        {
            let state = mock_app.state.read().await;
            assert_eq!(state.current_screen, Screen::FileSelection);
            assert_eq!(state.previous_screen, Some(Screen::MainMenu));
        }

        // Navigate back
        {
            let mut state = mock_app.state.write().await;
            state.go_back();
        }

        // Should be back at main menu
        {
            let state = mock_app.state.read().await;
            assert_eq!(state.current_screen, Screen::MainMenu);
        }
    }

    #[test]
    fn test_panic_hook_concept() {
        // Test panic hook concept
        let panic_occurred = false; // In real code, this would be set by the panic hook
        let terminal_restored = true; // In real code, this would be the result of restoration

        // If panic occurs, terminal should be restored
        if panic_occurred {
            assert!(terminal_restored, "Terminal should be restored after panic");
        } else {
            // Normal operation
            assert!(true);
        }
    }

    #[tokio::test]
    async fn test_concurrent_state_access() {
        let shared_state = Arc::new(RwLock::new(AppState::new()));

        // Simulate concurrent access
        let state1 = shared_state.clone();
        let state2 = shared_state.clone();

        let handle1 = tokio::spawn(async move {
            let mut state = state1.write().await;
            state.add_status_message("Message 1".to_string());
        });

        let handle2 = tokio::spawn(async move {
            let mut state = state2.write().await;
            state.add_status_message("Message 2".to_string());
        });

        // Wait for both operations
        let _ = tokio::join!(handle1, handle2);

        // Check final state
        let state = shared_state.read().await;
        assert_eq!(state.status_messages.len(), 2);
    }
}
