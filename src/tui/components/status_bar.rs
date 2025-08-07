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