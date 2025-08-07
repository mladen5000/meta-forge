use anyhow::Result;
use crossterm::event::{KeyCode, KeyEvent};
use ratatui::{
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    widgets::{Block, Borders, Paragraph},
    Frame,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use crate::tui::state::{AppState, Screen as ScreenEnum};
use super::Screen;

pub struct ErrorScreen {
    state: Arc<RwLock<AppState>>,
}

impl ErrorScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state }
    }
}

impl Screen for ErrorScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),  // Title
                Constraint::Min(5),     // Error message
                Constraint::Length(3),  // Actions
                Constraint::Length(2),  // Help
            ])
            .split(area);

        // Title
        let title = Paragraph::new("Error")
            .style(Style::default().fg(Color::Red).add_modifier(Modifier::BOLD))
            .alignment(Alignment::Center)
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        // Error message - placeholder
        let error_text = "An error occurred. Please try again or return to the main menu.".to_string();

        let error_msg = Paragraph::new(error_text)
            .style(Style::default().fg(Color::Red))
            .block(Block::default().title("Error Details").borders(Borders::ALL))
            .wrap(ratatui::widgets::Wrap { trim: true });

        f.render_widget(error_msg, chunks[1]);

        // Action suggestions
        let actions = Paragraph::new("• Press Enter to acknowledge and return\n• Press R to retry the last operation\n• Press Esc to return to main menu")
            .style(Style::default().fg(Color::Yellow))
            .block(Block::default().title("Available Actions").borders(Borders::ALL));

        f.render_widget(actions, chunks[2]);

        // Help
        let help = Paragraph::new("Enter: Acknowledge, R: Retry, Esc: Back to menu")
            .style(Style::default().fg(Color::Gray))
            .alignment(Alignment::Center);

        f.render_widget(help, chunks[3]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Enter | KeyCode::Esc => {
                // Return to main menu
                Ok(Some(ScreenEnum::MainMenu))
            }
            KeyCode::Char('r') | KeyCode::Char('R') => {
                // Return to main menu for retry
                Ok(Some(ScreenEnum::MainMenu))
            }
            _ => Ok(None),
        }
    }
}