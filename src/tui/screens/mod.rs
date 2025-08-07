pub mod main_menu;
pub mod file_selection;
pub mod configuration;
pub mod analysis;
pub mod database;
pub mod results;
pub mod help;
pub mod error;

use anyhow::Result;
use crossterm::event::KeyEvent;
use ratatui::Frame;
use std::sync::Arc;
use tokio::sync::RwLock;

use crate::tui::state::{AppState, Screen as ScreenEnum};

/// Trait for all screen implementations
pub trait Screen {
    /// Render the screen
    fn render(&mut self, f: &mut Frame, area: ratatui::layout::Rect) -> Result<()>;
    
    /// Handle key input
    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>>;
    
    /// Update the screen on tick
    fn tick(&mut self) -> Result<()> {
        Ok(())
    }
    
    /// Initialize the screen
    fn initialize(&mut self) -> Result<()> {
        Ok(())
    }
}

/// Screen manager for handling different screen types
pub struct ScreenManager {
    state: Arc<RwLock<AppState>>,
    main_menu: main_menu::MainMenuScreen,
    file_selection: file_selection::FileSelectionScreen,
    configuration: configuration::ConfigurationScreen,
    analysis: analysis::AnalysisScreen,
    database: database::DatabaseScreen,
    results: results::ResultsScreen,
    help: help::HelpScreen,
    error: error::ErrorScreen,
}

impl ScreenManager {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self {
            state: state.clone(),
            main_menu: main_menu::MainMenuScreen::new(state.clone()),
            file_selection: file_selection::FileSelectionScreen::new(state.clone()),
            configuration: configuration::ConfigurationScreen::new(state.clone()),
            analysis: analysis::AnalysisScreen::new(state.clone()),
            database: database::DatabaseScreen::new(state.clone()),
            results: results::ResultsScreen::new(state.clone()),
            help: help::HelpScreen::new(state.clone()),
            error: error::ErrorScreen::new(state.clone()),
        }
    }
    
    pub async fn initialize(&mut self) -> Result<()> {
        self.main_menu.initialize()?;
        self.file_selection.initialize()?;
        self.configuration.initialize()?;
        self.analysis.initialize()?;
        self.database.initialize()?;
        self.results.initialize()?;
        self.help.initialize()?;
        self.error.initialize()?;
        Ok(())
    }
    
    pub fn render(&mut self, f: &mut Frame, current_screen: &ScreenEnum) -> Result<()> {
        let area = f.area();
        
        match current_screen {
            ScreenEnum::MainMenu => self.main_menu.render(f, area),
            ScreenEnum::FileSelection => self.file_selection.render(f, area),
            ScreenEnum::Configuration => self.configuration.render(f, area),
            ScreenEnum::Analysis => self.analysis.render(f, area),
            ScreenEnum::Database => self.database.render(f, area),
            ScreenEnum::Results => self.results.render(f, area),
            ScreenEnum::Help => self.help.render(f, area),
            ScreenEnum::Error(_) => self.error.render(f, area),
        }
    }
    
    pub async fn handle_key_input(&mut self, key: KeyEvent, current_screen: &ScreenEnum) -> Result<()> {
        let new_screen = match current_screen {
            ScreenEnum::MainMenu => self.main_menu.handle_key(key)?,
            ScreenEnum::FileSelection => self.file_selection.handle_key(key)?,
            ScreenEnum::Configuration => self.configuration.handle_key(key)?,
            ScreenEnum::Analysis => self.analysis.handle_key(key)?,
            ScreenEnum::Database => self.database.handle_key(key)?,
            ScreenEnum::Results => self.results.handle_key(key)?,
            ScreenEnum::Help => self.help.handle_key(key)?,
            ScreenEnum::Error(_) => self.error.handle_key(key)?,
        };
        
        if let Some(screen) = new_screen {
            let mut state = self.state.write().await;
            state.navigate_to(screen);
        }
        
        Ok(())
    }
    
    pub async fn tick(&mut self) -> Result<()> {
        self.main_menu.tick()?;
        self.file_selection.tick()?;
        self.configuration.tick()?;
        self.analysis.tick()?;
        self.database.tick()?;
        self.results.tick()?;
        self.help.tick()?;
        self.error.tick()?;
        Ok(())
    }
}