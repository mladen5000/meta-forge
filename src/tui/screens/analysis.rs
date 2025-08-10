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

use super::Screen;
use crate::tui::state::{AppState, Screen as ScreenEnum};

pub struct AnalysisScreen {
    state: Arc<RwLock<AppState>>,
}

impl AnalysisScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state }
    }
}

impl Screen for AnalysisScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Title
                Constraint::Min(8),    // Progress
                Constraint::Length(3), // Status
                Constraint::Length(2), // Help
            ])
            .split(area);

        let title = Paragraph::new("Analysis Progress")
            .style(
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            )
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        // Progress display - placeholder for now
        let no_analysis = Paragraph::new(
            "No analysis running.\n\nUse 'Start Analysis' to begin processing your data.",
        )
        .style(Style::default().fg(Color::Yellow))
        .block(Block::default().title("Status").borders(Borders::ALL));
        f.render_widget(no_analysis, chunks[1]);

        // Status
        let status_text = "Ready to start analysis";

        let status = Paragraph::new(status_text)
            .style(Style::default().fg(Color::Green))
            .block(Block::default().title("Status").borders(Borders::ALL));
        f.render_widget(status, chunks[2]);

        let help = Paragraph::new("Space: Pause/Resume, Ctrl+C: Cancel, R: Restart, Esc: Back");
        f.render_widget(help, chunks[3]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Enter => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Char('r') | KeyCode::Char('R') => {
                // Restart analysis - placeholder
                Ok(None)
            }
            _ => Ok(None),
        }
    }
}
