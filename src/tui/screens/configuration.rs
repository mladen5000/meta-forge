use anyhow::Result;
use crossterm::event::{KeyCode, KeyEvent};
use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    widgets::{Block, Borders, Paragraph},
    Frame,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use crate::tui::state::{AppState, Screen as ScreenEnum};
use super::Screen;

pub struct ConfigurationScreen {
    state: Arc<RwLock<AppState>>,
}

impl ConfigurationScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state }
    }
}

impl Screen for ConfigurationScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),  // Title
                Constraint::Min(5),     // Config options
                Constraint::Length(2),  // Help
            ])
            .split(area);

        let title = Paragraph::new("Analysis Configuration")
            .style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        let config_text = "Sample Name: [sample_1]\nK-mer Range: [21-127]\nMin Coverage: [2]\nThreads: [auto]\nMemory Limit: [8GB]";
        let config = Paragraph::new(config_text)
            .block(Block::default().title("Configuration").borders(Borders::ALL));
        f.render_widget(config, chunks[1]);

        let help = Paragraph::new("↑/↓: Navigate, Enter: Edit, Tab: Next field, Esc: Back");
        f.render_widget(help, chunks[2]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Enter => Ok(Some(ScreenEnum::MainMenu)),
            _ => Ok(None),
        }
    }
}