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
                Constraint::Length(3),  // Title
                Constraint::Min(5),     // Results
                Constraint::Length(2),  // Help
            ])
            .split(area);

        let title = Paragraph::new("Analysis Results")
            .style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        let rt = tokio::runtime::Handle::current();
        let state = rt.block_on(async { self.state.read().await });
        
        let results_text = if state.results.sample_name.is_empty() {
            "No results available.\n\nRun an analysis to view results here.".to_string()
        } else {
            format!(
                "Sample: {}\nContigs: {}\nTotal Length: {} bp\nN50: {} bp\nProcessing Time: {:.2}s",
                state.results.sample_name,
                state.results.contigs_count,
                state.results.total_length,
                state.results.n50,
                state.results.processing_time.as_secs_f64()
            )
        };

        let results = Paragraph::new(results_text)
            .block(Block::default().title("Results Summary").borders(Borders::ALL));
        f.render_widget(results, chunks[1]);

        let help = Paragraph::new("E: Export, S: Save report, V: Visualize, Esc: Back");
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