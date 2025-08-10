use anyhow::Result;
use crossterm::event::{KeyCode, KeyEvent};
use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    widgets::{Block, Borders, List, ListItem, Paragraph},
    Frame,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use super::Screen;
use crate::tui::state::{AppState, Screen as ScreenEnum};

pub struct FileSelectionScreen {
    state: Arc<RwLock<AppState>>,
    selected_index: usize,
}

impl FileSelectionScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self {
            state,
            selected_index: 0,
        }
    }
}

impl Screen for FileSelectionScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Title
                Constraint::Min(5),    // File list
                Constraint::Length(3), // Selected files
                Constraint::Length(2), // Help
            ])
            .split(area);

        // Title
        let title = Paragraph::new("File Selection")
            .style(
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            )
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        // File list (placeholder)
        let files = vec![
            ListItem::new("sample_1.fastq"),
            ListItem::new("sample_2.fastq"),
            ListItem::new("reads.fasta"),
        ];

        let file_list = List::new(files)
            .block(
                Block::default()
                    .title("Available Files")
                    .borders(Borders::ALL),
            )
            .highlight_style(Style::default().bg(Color::DarkGray));

        f.render_widget(file_list, chunks[1]);

        // Selected files
        let selected = Paragraph::new("No files selected").block(
            Block::default()
                .title("Selected Files")
                .borders(Borders::ALL),
        );
        f.render_widget(selected, chunks[2]);

        // Help
        let help = Paragraph::new("↑/↓: Navigate, Space: Select, Enter: Confirm, Esc: Back");
        f.render_widget(help, chunks[3]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Up => {
                if self.selected_index > 0 {
                    self.selected_index -= 1;
                }
                Ok(None)
            }
            KeyCode::Down => {
                // Max 3 files in placeholder
                if self.selected_index < 2 {
                    self.selected_index += 1;
                }
                Ok(None)
            }
            KeyCode::Enter => Ok(Some(ScreenEnum::MainMenu)),
            _ => Ok(None),
        }
    }
}
