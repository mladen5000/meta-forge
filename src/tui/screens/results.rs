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

pub struct ResultsScreen {
    state: Arc<RwLock<AppState>>,
}

impl ResultsScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state }
    }
}

impl Screen for ResultsScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Title
                Constraint::Min(5),    // Results
                Constraint::Length(2), // Help
            ])
            .split(area);

        let title = Paragraph::new("Analysis Results")
            .style(
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            )
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        let results_text =
            "No results available.\n\nRun an analysis to view results here.".to_string();

        let results = Paragraph::new(results_text).block(
            Block::default()
                .title("Results Summary")
                .borders(Borders::ALL),
        );
        f.render_widget(results, chunks[1]);

        let help = Paragraph::new("E: Export, S: Save report, V: Visualize, Esc: Back");
        f.render_widget(help, chunks[2]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Enter => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Char('e') | KeyCode::Char('E') => {
                // Export results - placeholder
                Ok(None)
            }
            KeyCode::Char('s') | KeyCode::Char('S') => {
                // Save report - placeholder
                Ok(None)
            }
            KeyCode::Char('v') | KeyCode::Char('V') => {
                // Visualize data - placeholder
                Ok(None)
            }
            _ => Ok(None),
        }
    }
}
