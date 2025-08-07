use anyhow::Result;
use crossterm::event::{KeyCode, KeyEvent};
use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, List, ListItem, Paragraph},
    Frame,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use crate::tui::state::{AppState, Screen as ScreenEnum};
use super::Screen;

pub struct HelpScreen {
    state: Arc<RwLock<AppState>>,
    scroll: usize,
}

impl HelpScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        Self { state, scroll: 0 }
    }
}

impl Screen for HelpScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),  // Title
                Constraint::Min(5),     // Help content
                Constraint::Length(2),  // Navigation help
            ])
            .split(area);

        let title = Paragraph::new("Help & Documentation")
            .style(Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        // Help sections
        let help_items = vec![
            ListItem::new(Line::from(vec![
                Span::styled("General Navigation:", Style::default().add_modifier(Modifier::BOLD)),
            ])),
            ListItem::new("  ↑/↓ or j/k - Navigate up/down"),
            ListItem::new("  Enter - Select/confirm"),
            ListItem::new("  Esc - Go back/cancel"),
            ListItem::new("  Ctrl+Q - Quit application"),
            ListItem::new("  F1 - Show this help"),
            ListItem::new(""),
            ListItem::new(Line::from(vec![
                Span::styled("Main Menu:", Style::default().add_modifier(Modifier::BOLD)),
            ])),
            ListItem::new("  1-6 - Quick selection by number"),
            ListItem::new("  h - Show help"),
            ListItem::new("  q - Quit"),
            ListItem::new(""),
            ListItem::new(Line::from(vec![
                Span::styled("File Selection:", Style::default().add_modifier(Modifier::BOLD)),
            ])),
            ListItem::new("  Space - Select/deselect file"),
            ListItem::new("  Ctrl+A - Select all files"),
            ListItem::new("  Ctrl+D - Deselect all files"),
            ListItem::new("  F5 - Refresh file list"),
            ListItem::new(""),
            ListItem::new(Line::from(vec![
                Span::styled("Analysis:", Style::default().add_modifier(Modifier::BOLD)),
            ])),
            ListItem::new("  Space - Pause/resume operation"),
            ListItem::new("  Ctrl+C - Cancel current operation"),
            ListItem::new("  R - Restart analysis"),
            ListItem::new(""),
            ListItem::new(Line::from(vec![
                Span::styled("About Meta-Forge:", Style::default().add_modifier(Modifier::BOLD)),
            ])),
            ListItem::new("  A comprehensive metagenomic analysis toolkit"),
            ListItem::new("  with integrated machine learning capabilities"),
            ListItem::new("  Version: 0.1.0"),
        ];

        let help_list = List::new(help_items)
            .block(Block::default().title("Help Topics").borders(Borders::ALL));

        f.render_widget(help_list, chunks[1]);

        let nav_help = Paragraph::new("↑/↓: Scroll, Page Up/Down: Page scroll, Esc: Back to menu");
        f.render_widget(nav_help, chunks[2]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Esc => Ok(Some(ScreenEnum::MainMenu)),
            KeyCode::Up => {
                if self.scroll > 0 {
                    self.scroll -= 1;
                }
                Ok(None)
            }
            KeyCode::Down => {
                self.scroll += 1;
                Ok(None)
            }
            KeyCode::PageUp => {
                self.scroll = self.scroll.saturating_sub(10);
                Ok(None)
            }
            KeyCode::PageDown => {
                self.scroll += 10;
                Ok(None)
            }
            _ => Ok(None),
        }
    }
}