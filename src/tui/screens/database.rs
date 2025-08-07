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

pub struct DatabaseScreen {
    state: Arc<RwLock<AppState>>,
}

impl DatabaseScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state }
    }
}

impl Screen for DatabaseScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),  // Title
                Constraint::Min(5),     // Database info
                Constraint::Length(2),  // Help
            ])
            .split(area);

        let title = Paragraph::new("Database Management")
            .style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        let rt = tokio::runtime::Handle::current();
        let state = rt.block_on(async { self.state.read().await });
        
        let db_info = if state.database_info.initialized {
            format!(
                "Database Status: Initialized\nPath: {:?}\nSize: {:?}\nTables: {:?}",
                state.database_info.path,
                state.database_info.size,
                state.database_info.tables.len()
            )
        } else {
            "Database Status: Not initialized\n\nPress '1' to initialize the database".to_string()
        };

        let db_panel = Paragraph::new(db_info)
            .block(Block::default().title("Database Information").borders(Borders::ALL));
        f.render_widget(db_panel, chunks[1]);

        let help = Paragraph::new("1: Initialize, 2: Backup, 3: Restore, Esc: Back");
        f.render_widget(help, chunks[2]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            _ => Ok(None),
        }
    }
}