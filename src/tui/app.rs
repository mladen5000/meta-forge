use anyhow::{Context, Result};
use crossterm::{
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::CrosstermBackend,
    Terminal,
};
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
                eprintln!("Render error: {}", e);
            }
        })?;
        
        Ok(())
    }
    
    /// Handle application events
    async fn handle_event(&mut self, event: AppEvent) -> Result<()> {
        match event {
            AppEvent::Input(key) => {
                // Handle global key bindings first
                match key.code {
                    crossterm::event::KeyCode::Char('q') => {
                        if key.modifiers.contains(crossterm::event::KeyModifiers::CONTROL) {
                            let mut state = self.state.write().await;
                            state.should_quit = true;
                            return Ok(());
                        }
                    }
                    crossterm::event::KeyCode::Esc => {
                        let state = self.state.read().await;
                        if matches!(state.current_screen, crate::tui::state::Screen::MainMenu) {
                            drop(state);
                            let mut state = self.state.write().await;
                            state.should_quit = true;
                            return Ok(());
                        }
                        drop(state);
                    }
                    _ => {}
                }
                
                // Get current screen for handling
                let screen = {
                    let state = self.state.read().await;
                    state.current_screen.clone()
                };
                
                self.screen_manager.handle_key_input(key, &screen).await?;
            }
            
            AppEvent::Tick => {
                // Update timers and refresh displays
                self.screen_manager.tick().await?;
            }
            
            AppEvent::AnalysisStarted(name) => {
                let mut state = self.state.write().await;
                state.add_operation(&name, None);
                state.add_status_message(format!("Started: {}", name));
            }
            
            AppEvent::AnalysisProgress { operation, current, total, message } => {
                let mut state = self.state.write().await;
                if !state.operations.contains_key(&operation) {
                    state.add_operation(&operation, total);
                }
                state.update_operation(&operation, current, Some(&message));
            }
            
            AppEvent::AnalysisCompleted(name) => {
                let mut state = self.state.write().await;
                state.complete_operation(&name);
                state.add_status_message(format!("Completed: {}", name));
            }
            
            AppEvent::AnalysisError(name, error) => {
                let mut state = self.state.write().await;
                state.error_operation(&name, &error);
                state.set_error(format!("Analysis error in {}: {}", name, error));
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
        let _ = execute!(
            self.terminal.backend_mut(),
            LeaveAlternateScreen
        );
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
        eprintln!("Application panicked: {}", panic);
    }));
    
    let result = app.run().await;
    
    // Restore default panic hook
    let _ = std::panic::take_hook();
    
    result
}