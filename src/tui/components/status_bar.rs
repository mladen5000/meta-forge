use ratatui::{
    layout::Rect,
    style::{Color, Style},
    widgets::{Block, Borders, Paragraph},
    Frame,
};

/// Status bar component for displaying system status
pub struct StatusBar {
    message: String,
    style: Style,
}

impl StatusBar {
    pub fn new(message: String) -> Self {
        Self {
            message,
            style: Style::default().fg(Color::Green),
        }
    }

    pub fn with_style(mut self, style: Style) -> Self {
        self.style = style;
        self
    }

    pub fn error(message: String) -> Self {
        Self {
            message,
            style: Style::default().fg(Color::Red),
        }
    }

    pub fn warning(message: String) -> Self {
        Self {
            message,
            style: Style::default().fg(Color::Yellow),
        }
    }

    pub fn render(&self, f: &mut Frame, area: Rect) {
        let status = Paragraph::new(self.message.as_str())
            .style(self.style)
            .block(Block::default().borders(Borders::ALL));

        f.render_widget(status, area);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ratatui::style::Color;

    #[test]
    fn test_status_bar_new() {
        let status_bar = StatusBar::new("Test message".to_string());

        assert_eq!(status_bar.message, "Test message");
        assert_eq!(status_bar.style, Style::default().fg(Color::Green));
    }

    #[test]
    fn test_status_bar_with_style() {
        let custom_style = Style::default().fg(Color::Blue);
        let status_bar = StatusBar::new("Test".to_string()).with_style(custom_style);

        assert_eq!(status_bar.message, "Test");
        assert_eq!(status_bar.style, custom_style);
    }

    #[test]
    fn test_status_bar_error() {
        let error_bar = StatusBar::error("Error occurred".to_string());

        assert_eq!(error_bar.message, "Error occurred");
        assert_eq!(error_bar.style, Style::default().fg(Color::Red));
    }

    #[test]
    fn test_status_bar_warning() {
        let warning_bar = StatusBar::warning("Warning message".to_string());

        assert_eq!(warning_bar.message, "Warning message");
        assert_eq!(warning_bar.style, Style::default().fg(Color::Yellow));
    }

    #[test]
    fn test_status_bar_style_chaining() {
        let status_bar =
            StatusBar::new("Test".to_string()).with_style(Style::default().fg(Color::Cyan));

        assert_eq!(status_bar.style.fg, Some(Color::Cyan));
    }

    #[test]
    fn test_empty_message() {
        let status_bar = StatusBar::new(String::new());
        assert!(status_bar.message.is_empty());
    }

    #[test]
    fn test_long_message() {
        let long_message = "A".repeat(1000);
        let status_bar = StatusBar::new(long_message.clone());
        assert_eq!(status_bar.message, long_message);
    }

    #[test]
    fn test_unicode_message() {
        let unicode_message = "🧬 DNA Analysis Complete ✅".to_string();
        let status_bar = StatusBar::new(unicode_message.clone());
        assert_eq!(status_bar.message, unicode_message);
    }

    #[test]
    fn test_multiline_message() {
        let multiline_message = "Line 1\nLine 2\nLine 3".to_string();
        let status_bar = StatusBar::new(multiline_message.clone());
        assert_eq!(status_bar.message, multiline_message);
    }

    #[test]
    fn test_status_bar_methods_return_correct_colors() {
        let normal = StatusBar::new("normal".to_string());
        let error = StatusBar::error("error".to_string());
        let warning = StatusBar::warning("warning".to_string());

        assert_eq!(normal.style.fg, Some(Color::Green));
        assert_eq!(error.style.fg, Some(Color::Red));
        assert_eq!(warning.style.fg, Some(Color::Yellow));
    }
}
