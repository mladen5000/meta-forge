use anyhow::Result;
use crossterm::event::{self, Event, KeyCode, KeyEvent, KeyModifiers};
use std::sync::mpsc;
use std::thread;
use std::time::{Duration, Instant};
use tokio::sync::mpsc as tokio_mpsc;

/// Application events
#[derive(Debug, Clone)]
pub enum AppEvent {
    /// Terminal input event
    Input(KeyEvent),
    
    /// Timer tick for updates
    Tick,
    
    /// Analysis started
    AnalysisStarted(String),
    
    /// Analysis progress update
    AnalysisProgress {
        operation: String,
        current: u64,
        total: Option<u64>,
        message: String,
    },
    
    /// Analysis completed
    AnalysisCompleted(String),
    
    /// Analysis failed
    AnalysisError(String, String),
    
    /// Database status update
    DatabaseUpdate(crate::tui::state::DatabaseInfo),
    
    /// File list updated
    FileListUpdated(Vec<std::path::PathBuf>),
    
    /// Results available
    ResultsAvailable(crate::tui::state::ResultsData),
    
    /// Status message
    StatusMessage(String),
    
    /// Quit application
    Quit,
}

/// Event handler for the TUI application
pub struct EventHandler {
    /// Receiver for input events
    receiver: mpsc::Receiver<AppEvent>,
    
    /// Sender for input events
    _sender: mpsc::Sender<AppEvent>,
    
    /// Async event sender
    async_sender: tokio_mpsc::UnboundedSender<AppEvent>,
    
    /// Async event receiver
    pub async_receiver: tokio_mpsc::UnboundedReceiver<AppEvent>,
    
    /// Handler thread handle
    _handle: thread::JoinHandle<()>,
}

impl EventHandler {
    pub fn new(tick_rate: Duration) -> Self {
        let (sender, receiver) = mpsc::channel();
        let (async_sender, async_receiver) = tokio_mpsc::unbounded_channel();
        
        let event_sender = sender.clone();
        let handle = thread::spawn(move || {
            let mut last_tick = Instant::now();
            
            loop {
                let timeout = tick_rate
                    .checked_sub(last_tick.elapsed())
                    .unwrap_or_else(|| Duration::from_secs(0));
                
                if event::poll(timeout).unwrap_or(false) {
                    match event::read() {
                        Ok(Event::Key(key)) => {
                            if event_sender.send(AppEvent::Input(key)).is_err() {
                                break;
                            }
                        }
                        Ok(_) => {}
                        Err(_) => break,
                    }
                }
                
                if last_tick.elapsed() >= tick_rate {
                    if event_sender.send(AppEvent::Tick).is_err() {
                        break;
                    }
                    last_tick = Instant::now();
                }
            }
        });
        
        Self {
            receiver,
            _sender: sender,
            async_sender,
            async_receiver,
            _handle: handle,
        }
    }
    
    /// Get the next event
    pub fn next(&self) -> Result<AppEvent> {
        Ok(self.receiver.recv()?)
    }
    
    /// Try to get an event without blocking
    pub fn try_next(&self) -> Result<Option<AppEvent>> {
        match self.receiver.try_recv() {
            Ok(event) => Ok(Some(event)),
            Err(mpsc::TryRecvError::Empty) => Ok(None),
            Err(mpsc::TryRecvError::Disconnected) => {
                Err(anyhow::anyhow!("Event channel disconnected"))
            }
        }
    }
    
    /// Get async event sender for background tasks
    pub fn async_sender(&self) -> tokio_mpsc::UnboundedSender<AppEvent> {
        self.async_sender.clone()
    }
}

/// Handle key events and return appropriate application events
pub fn handle_key_event(key: KeyEvent, current_screen: &crate::tui::state::Screen) -> Option<AppEvent> {
    use crate::tui::state::Screen;
    use KeyCode::*;
    
    match key.code {
        // Global key bindings
        Char('q') | Char('Q') => {
            if key.modifiers.contains(KeyModifiers::CONTROL) {
                Some(AppEvent::Quit)
            } else {
                None
            }
        }
        
        // Escape - go back or quit from main menu
        Esc => match current_screen {
            Screen::MainMenu => Some(AppEvent::Quit),
            _ => None, // Let the screen handle going back
        },
        
        // Screen-specific navigation
        Char('h') | F(1) => match current_screen {
            Screen::Help => None,
            _ => None, // Navigate to help - handled by screens
        },
        
        _ => None,
    }
}

/// Key binding information for help screen
pub struct KeyBinding {
    pub key: String,
    pub description: String,
}

impl KeyBinding {
    pub fn new(key: &str, description: &str) -> Self {
        Self {
            key: key.to_string(),
            description: description.to_string(),
        }
    }
}

/// Get key bindings for the current screen
pub fn get_key_bindings(screen: &crate::tui::state::Screen) -> Vec<KeyBinding> {
    use crate::tui::state::Screen;
    
    let mut bindings = vec![
        KeyBinding::new("Ctrl+Q", "Quit application"),
        KeyBinding::new("Esc", "Go back / Exit"),
        KeyBinding::new("F1", "Help"),
    ];
    
    match screen {
        Screen::MainMenu => {
            bindings.extend(vec![
                KeyBinding::new("1", "File Selection"),
                KeyBinding::new("2", "Configuration"),
                KeyBinding::new("3", "Start Analysis"),
                KeyBinding::new("4", "Database Management"),
                KeyBinding::new("5", "View Results"),
                KeyBinding::new("↑/↓", "Navigate menu"),
                KeyBinding::new("Enter", "Select item"),
            ]);
        }
        
        Screen::FileSelection => {
            bindings.extend(vec![
                KeyBinding::new("↑/↓", "Navigate files"),
                KeyBinding::new("Enter", "Select/Deselect file"),
                KeyBinding::new("Space", "Preview file"),
                KeyBinding::new("Tab", "Toggle filter"),
                KeyBinding::new("F5", "Refresh"),
                KeyBinding::new("Ctrl+A", "Select all"),
                KeyBinding::new("Ctrl+D", "Deselect all"),
            ]);
        }
        
        Screen::Configuration => {
            bindings.extend(vec![
                KeyBinding::new("↑/↓", "Navigate fields"),
                KeyBinding::new("Enter", "Edit field"),
                KeyBinding::new("Tab", "Next field"),
                KeyBinding::new("Shift+Tab", "Previous field"),
                KeyBinding::new("F2", "Save configuration"),
                KeyBinding::new("F3", "Load defaults"),
            ]);
        }
        
        Screen::Analysis => {
            bindings.extend(vec![
                KeyBinding::new("↑/↓", "Navigate operations"),
                KeyBinding::new("Space", "Pause/Resume"),
                KeyBinding::new("Ctrl+C", "Cancel operation"),
                KeyBinding::new("F5", "Refresh status"),
                KeyBinding::new("Enter", "View details"),
            ]);
        }
        
        Screen::Database => {
            bindings.extend(vec![
                KeyBinding::new("1", "Initialize database"),
                KeyBinding::new("2", "Check status"),
                KeyBinding::new("3", "Backup database"),
                KeyBinding::new("4", "Restore database"),
                KeyBinding::new("F5", "Refresh status"),
            ]);
        }
        
        Screen::Results => {
            bindings.extend(vec![
                KeyBinding::new("↑/↓", "Navigate results"),
                KeyBinding::new("Enter", "View details"),
                KeyBinding::new("E", "Export results"),
                KeyBinding::new("V", "Visualize data"),
                KeyBinding::new("S", "Save report"),
            ]);
        }
        
        Screen::Help => {
            bindings.extend(vec![
                KeyBinding::new("↑/↓", "Scroll help"),
                KeyBinding::new("Page Up/Down", "Page scroll"),
            ]);
        }
        
        Screen::Error(_) => {
            bindings.extend(vec![
                KeyBinding::new("Enter", "Acknowledge error"),
                KeyBinding::new("R", "Retry operation"),
            ]);
        }
    }
    
    bindings
}