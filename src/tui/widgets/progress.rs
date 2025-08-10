use ratatui::{
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Gauge, List, ListItem, Paragraph},
    Frame,
};
use std::collections::HashMap;

use crate::tui::state::{OperationProgress, OperationStatus};

/// Progress bar widget for single operations
pub struct ProgressWidget<'a> {
    progress: &'a OperationProgress,
    title: Option<&'a str>,
}

impl<'a> ProgressWidget<'a> {
    pub fn new(progress: &'a OperationProgress) -> Self {
        Self {
            progress,
            title: None,
        }
    }

    pub fn title(mut self, title: &'a str) -> Self {
        self.title = Some(title);
        self
    }

    pub fn render(&self, f: &mut Frame, area: Rect) {
        let block = Block::default()
            .title(self.title.unwrap_or(&self.progress.name))
            .borders(Borders::ALL);

        let inner = block.inner(area);
        f.render_widget(block, area);

        // Split area for progress bar and info
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Progress bar
                Constraint::Min(1),    // Info
            ])
            .split(inner);

        // Render progress bar
        let percentage = self.progress.percentage().unwrap_or(0.0);
        let gauge = Gauge::default()
            .block(Block::default().borders(Borders::NONE))
            .gauge_style(match self.progress.status {
                OperationStatus::Completed => Style::default().fg(Color::Green),
                OperationStatus::Error(_) => Style::default().fg(Color::Red),
                OperationStatus::Running => Style::default().fg(Color::Blue),
                OperationStatus::Cancelled => Style::default().fg(Color::Yellow),
                OperationStatus::Pending => Style::default().fg(Color::Gray),
            })
            .percent(percentage as u16)
            .label(format!("{percentage:.1}%"));

        f.render_widget(gauge, chunks[0]);

        // Render info
        let mut lines = vec![];

        // Status line
        let status_text = match &self.progress.status {
            OperationStatus::Pending => "Pending",
            OperationStatus::Running => "Running",
            OperationStatus::Completed => "Completed",
            OperationStatus::Error(e) => &format!("Error: {e}"),
            OperationStatus::Cancelled => "Cancelled",
        };

        lines.push(Line::from(vec![
            Span::styled("Status: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::styled(
                status_text,
                Style::default().fg(match self.progress.status {
                    OperationStatus::Completed => Color::Green,
                    OperationStatus::Error(_) => Color::Red,
                    OperationStatus::Running => Color::Blue,
                    OperationStatus::Cancelled => Color::Yellow,
                    OperationStatus::Pending => Color::Gray,
                }),
            ),
        ]));

        // Progress counts
        if let Some(total) = self.progress.total {
            lines.push(Line::from(vec![
                Span::styled("Progress: ", Style::default().add_modifier(Modifier::BOLD)),
                Span::raw(format!("{} / {} items", self.progress.current, total)),
            ]));
        } else {
            lines.push(Line::from(vec![
                Span::styled("Progress: ", Style::default().add_modifier(Modifier::BOLD)),
                Span::raw(format!("{} items", self.progress.current)),
            ]));
        }

        // Rate and ETA
        if let Some(rate) = self.progress.rate {
            let rate_text = if rate >= 1000.0 {
                format!("{:.1}K/s", rate / 1000.0)
            } else {
                format!("{rate:.0}/s")
            };

            let mut rate_line = vec![
                Span::styled("Rate: ", Style::default().add_modifier(Modifier::BOLD)),
                Span::raw(rate_text),
            ];

            if let Some(eta) = self.progress.eta {
                rate_line.extend(vec![Span::raw(" | ETA: "), Span::raw(format_duration(eta))]);
            }

            lines.push(Line::from(rate_line));
        }

        // Elapsed time
        let elapsed = self.progress.start_time.elapsed();
        lines.push(Line::from(vec![
            Span::styled("Elapsed: ", Style::default().add_modifier(Modifier::BOLD)),
            Span::raw(format_duration(elapsed)),
        ]));

        // Current message
        if !self.progress.message.is_empty() {
            lines.push(Line::from(vec![
                Span::styled("Message: ", Style::default().add_modifier(Modifier::BOLD)),
                Span::raw(&self.progress.message),
            ]));
        }

        let info = Paragraph::new(lines)
            .block(Block::default().borders(Borders::NONE))
            .wrap(ratatui::widgets::Wrap { trim: true });

        f.render_widget(info, chunks[1]);
    }
}

/// Multi-progress widget for multiple concurrent operations
pub struct MultiProgressWidget<'a> {
    operations: &'a HashMap<String, OperationProgress>,
    title: Option<&'a str>,
    selected: Option<usize>,
}

impl<'a> MultiProgressWidget<'a> {
    pub fn new(operations: &'a HashMap<String, OperationProgress>) -> Self {
        Self {
            operations,
            title: None,
            selected: None,
        }
    }

    pub fn title(mut self, title: &'a str) -> Self {
        self.title = Some(title);
        self
    }

    pub fn selected(mut self, selected: Option<usize>) -> Self {
        self.selected = selected;
        self
    }

    pub fn render(&self, f: &mut Frame, area: Rect) {
        let block = Block::default()
            .title(self.title.unwrap_or("Operations"))
            .borders(Borders::ALL);

        let inner = block.inner(area);
        f.render_widget(block, area);

        if self.operations.is_empty() {
            let no_ops = Paragraph::new("No active operations")
                .style(Style::default().fg(Color::Gray))
                .block(Block::default().borders(Borders::NONE));
            f.render_widget(no_ops, inner);
            return;
        }

        // Create list items for each operation
        let items: Vec<ListItem> = self
            .operations
            .iter()
            .enumerate()
            .map(|(i, (name, progress))| {
                let percentage = progress.percentage().unwrap_or(0.0);
                let status_symbol = match progress.status {
                    OperationStatus::Completed => "✓",
                    OperationStatus::Error(_) => "✗",
                    OperationStatus::Running => "⚡",
                    OperationStatus::Cancelled => "⚠",
                    OperationStatus::Pending => "⏸",
                };

                let status_color = match progress.status {
                    OperationStatus::Completed => Color::Green,
                    OperationStatus::Error(_) => Color::Red,
                    OperationStatus::Running => Color::Blue,
                    OperationStatus::Cancelled => Color::Yellow,
                    OperationStatus::Pending => Color::Gray,
                };

                let progress_bar = create_progress_bar(percentage, 20);

                let progress_text = if let Some(total) = progress.total {
                    format!("{:.1}% ({}/{})", percentage, progress.current, total)
                } else {
                    format!("{}", progress.current)
                };

                let style = if Some(i) == self.selected {
                    Style::default().bg(Color::DarkGray)
                } else {
                    Style::default()
                };

                ListItem::new(Line::from(vec![
                    Span::styled(status_symbol, Style::default().fg(status_color)),
                    Span::raw(" "),
                    Span::styled(name, style),
                    Span::raw(" ["),
                    Span::styled(progress_bar, Style::default().fg(status_color)),
                    Span::raw("] "),
                    Span::raw(progress_text),
                ]))
            })
            .collect();

        let list = List::new(items)
            .block(Block::default().borders(Borders::NONE))
            .highlight_style(Style::default().bg(Color::DarkGray))
            .highlight_symbol("► ");

        f.render_widget(list, inner);
    }
}

/// Compact operation progress widget for dashboard
pub struct OperationProgressWidget<'a> {
    progress: &'a OperationProgress,
}

impl<'a> OperationProgressWidget<'a> {
    pub fn new(progress: &'a OperationProgress) -> Self {
        Self { progress }
    }

    pub fn render(&self, f: &mut Frame, area: Rect) {
        let percentage = self.progress.percentage().unwrap_or(0.0);
        let progress_bar = create_progress_bar(percentage, area.width.saturating_sub(4) as usize);

        let status_color = match self.progress.status {
            OperationStatus::Completed => Color::Green,
            OperationStatus::Error(_) => Color::Red,
            OperationStatus::Running => Color::Blue,
            OperationStatus::Cancelled => Color::Yellow,
            OperationStatus::Pending => Color::Gray,
        };

        let text = if let Some(total) = self.progress.total {
            format!(
                "{} [{progress_bar}] {:.1}% ({}/{})",
                self.progress.name, percentage, self.progress.current, total
            )
        } else {
            format!(
                "{} [{progress_bar}] {}",
                self.progress.name, self.progress.current
            )
        };

        let paragraph = Paragraph::new(text)
            .style(Style::default().fg(status_color))
            .block(Block::default().borders(Borders::NONE));

        f.render_widget(paragraph, area);
    }
}

/// Create a text-based progress bar
fn create_progress_bar(percentage: f64, width: usize) -> String {
    let filled = ((percentage / 100.0) * width as f64) as usize;
    let empty = width.saturating_sub(filled);

    format!("{}{}", "█".repeat(filled), "░".repeat(empty))
}

/// Format duration in a human-readable way
fn format_duration(duration: std::time::Duration) -> String {
    let secs = duration.as_secs();

    if secs < 60 {
        format!("{secs}s")
    } else if secs < 3600 {
        format!("{}m {}s", secs / 60, secs % 60)
    } else {
        format!("{}h {}m", secs / 3600, (secs % 3600) / 60)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_progress_bar() {
        assert_eq!(create_progress_bar(0.0, 10), "░░░░░░░░░░");
        assert_eq!(create_progress_bar(50.0, 10), "█████░░░░░");
        assert_eq!(create_progress_bar(100.0, 10), "██████████");
    }

    #[test]
    fn test_format_duration() {
        use std::time::Duration;

        assert_eq!(format_duration(Duration::from_secs(30)), "30s");
        assert_eq!(format_duration(Duration::from_secs(90)), "1m 30s");
        assert_eq!(format_duration(Duration::from_secs(3661)), "1h 1m");
    }
}
