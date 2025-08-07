use ratatui::{
    layout::Rect,
    style::{Color, Style},
    widgets::{Block, Borders, Paragraph},
    Frame,
};

/// Navigation help component
pub struct Navigation {
    help_text: String,
}

impl Navigation {
    pub fn new(help_text: String) -> Self {
        Self { help_text }
    }
    
    pub fn render(&self, f: &mut Frame, area: Rect) {
        let nav = Paragraph::new(self.help_text.as_str())
            .style(Style::default().fg(Color::Gray))
            .block(Block::default().borders(Borders::NONE));
        
        f.render_widget(nav, area);
    }
}