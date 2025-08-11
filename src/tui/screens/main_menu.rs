use anyhow::Result;
use crossterm::event::{KeyCode, KeyEvent};
use ratatui::{
    layout::{Alignment, Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, List, ListItem, Paragraph},
    Frame,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use super::Screen;
use crate::tui::state::{AppState, Screen as ScreenEnum};

/// Main menu screen
pub struct MainMenuScreen {
    state: Arc<RwLock<AppState>>,
    selected_item: usize,
    menu_items: Vec<MenuItem>,
}

struct MenuItem {
    title: &'static str,
    description: &'static str,
    screen: ScreenEnum,
    enabled: bool,
}

impl MainMenuScreen {
    pub fn new(state: Arc<RwLock<AppState>>) -> Self {
        let menu_items = vec![
            MenuItem {
                title: "File Selection",
                description: "Select input files for analysis",
                screen: ScreenEnum::FileSelection,
                enabled: true,
            },
            MenuItem {
                title: "Configuration",
                description: "Configure analysis parameters",
                screen: ScreenEnum::Configuration,
                enabled: true,
            },
            MenuItem {
                title: "Run Analysis",
                description: "Start metagenomic analysis",
                screen: ScreenEnum::Analysis,
                enabled: true,
            },
            MenuItem {
                title: "Database",
                description: "Manage analysis database",
                screen: ScreenEnum::Database,
                enabled: true,
            },
            MenuItem {
                title: "Results",
                description: "View analysis results",
                screen: ScreenEnum::Results,
                enabled: true,
            },
            MenuItem {
                title: "Help",
                description: "View help and documentation",
                screen: ScreenEnum::Help,
                enabled: true,
            },
        ];

        Self {
            state,
            selected_item: 0,
            menu_items,
        }
    }
}

impl Screen for MainMenuScreen {
    fn render(&mut self, f: &mut Frame, area: Rect) -> Result<()> {
        // Main layout
        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3), // Title
                Constraint::Min(10),   // Menu
                Constraint::Length(3), // Status
                Constraint::Length(2), // Help
            ])
            .split(area);

        // Title
        let title = Paragraph::new("Meta-Forge - Metagenomic Analysis Suite")
            .style(
                Style::default()
                    .fg(Color::Cyan)
                    .add_modifier(Modifier::BOLD),
            )
            .alignment(Alignment::Center)
            .block(Block::default().borders(Borders::ALL));
        f.render_widget(title, chunks[0]);

        // Menu items
        let menu_layout = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([
                Constraint::Percentage(60), // Menu list
                Constraint::Percentage(40), // Description
            ])
            .split(chunks[1]);

        // Create menu items
        let items: Vec<ListItem> = self
            .menu_items
            .iter()
            .enumerate()
            .map(|(i, item)| {
                let style = if item.enabled {
                    if i == self.selected_item {
                        Style::default()
                            .fg(Color::Yellow)
                            .add_modifier(Modifier::BOLD)
                    } else {
                        Style::default().fg(Color::White)
                    }
                } else {
                    Style::default().fg(Color::DarkGray)
                };

                let prefix = if i == self.selected_item {
                    "► "
                } else {
                    "  "
                };
                ListItem::new(Line::from(vec![
                    Span::styled(prefix, style),
                    Span::styled(format!("{}. {}", i + 1, item.title), style),
                ]))
            })
            .collect();

        let menu_list = List::new(items)
            .block(Block::default().title("Main Menu").borders(Borders::ALL))
            .highlight_style(Style::default().bg(Color::DarkGray));

        f.render_widget(menu_list, menu_layout[0]);

        // Description panel
        let selected_item = &self.menu_items[self.selected_item];
        let description = Paragraph::new(vec![
            Line::from(vec![
                Span::styled("Selected: ", Style::default().add_modifier(Modifier::BOLD)),
                Span::styled(selected_item.title, Style::default().fg(Color::Cyan)),
            ]),
            Line::from(""),
            Line::from(selected_item.description),
            Line::from(""),
            if selected_item.enabled {
                Line::from(vec![
                    Span::styled("Status: ", Style::default().add_modifier(Modifier::BOLD)),
                    Span::styled("Available", Style::default().fg(Color::Green)),
                ])
            } else {
                Line::from(vec![
                    Span::styled("Status: ", Style::default().add_modifier(Modifier::BOLD)),
                    Span::styled("Disabled", Style::default().fg(Color::Red)),
                ])
            },
        ])
        .block(Block::default().title("Description").borders(Borders::ALL))
        .wrap(ratatui::widgets::Wrap { trim: true });

        f.render_widget(description, menu_layout[1]);

        // Status bar - simplified to avoid async in render
        let status_text = "Ready - Use arrow keys to navigate, Enter to select".to_string();

        let status = Paragraph::new(status_text)
            .style(Style::default().fg(Color::Green))
            .alignment(Alignment::Left)
            .block(Block::default().borders(Borders::ALL).title("Status"));

        f.render_widget(status, chunks[2]);

        // Help text
        let help_text = "Use ↑/↓ to navigate, Enter to select, Ctrl+Q to quit, F1 for help";
        let help = Paragraph::new(help_text)
            .style(Style::default().fg(Color::Gray))
            .alignment(Alignment::Center);

        f.render_widget(help, chunks[3]);

        Ok(())
    }

    fn handle_key(&mut self, key: KeyEvent) -> Result<Option<ScreenEnum>> {
        match key.code {
            KeyCode::Up => {
                if self.selected_item > 0 {
                    self.selected_item -= 1;
                } else {
                    self.selected_item = self.menu_items.len() - 1;
                }
            }

            KeyCode::Down => {
                if self.selected_item < self.menu_items.len() - 1 {
                    self.selected_item += 1;
                } else {
                    self.selected_item = 0;
                }
            }

            KeyCode::Enter => {
                let selected = &self.menu_items[self.selected_item];
                if selected.enabled {
                    return Ok(Some(selected.screen.clone()));
                }
            }

            KeyCode::Char(c) => {
                if let Some(digit) = c.to_digit(10) {
                    let index = (digit as usize).saturating_sub(1);
                    if index < self.menu_items.len() {
                        self.selected_item = index;
                        let selected = &self.menu_items[self.selected_item];
                        if selected.enabled {
                            return Ok(Some(selected.screen.clone()));
                        }
                    }
                }

                match c {
                    'q' | 'Q' => {
                        // Quit from main menu
                    }
                    'h' | 'H' => return Ok(Some(ScreenEnum::Help)),
                    _ => {}
                }
            }

            KeyCode::F(1) => return Ok(Some(ScreenEnum::Help)),
            KeyCode::Esc => {
                // Quit from main menu
            }

            _ => {}
        }

        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tui::state::AppState;
    use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
    use std::sync::Arc;
    use tokio::sync::RwLock;

    fn create_test_screen() -> MainMenuScreen {
        let state = Arc::new(RwLock::new(AppState::new()));
        MainMenuScreen::new(state)
    }

    #[test]
    fn test_main_menu_creation() {
        let screen = create_test_screen();

        assert_eq!(screen.selected_item, 0);
        assert_eq!(screen.menu_items.len(), 6); // File, Config, Analysis, DB, Results, Help
    }

    #[test]
    fn test_menu_items_structure() {
        let screen = create_test_screen();

        // Check that all expected menu items exist
        let expected_titles = vec![
            "File Selection",
            "Configuration",
            "Run Analysis",
            "Database",
            "Results",
            "Help",
        ];

        for (i, expected_title) in expected_titles.iter().enumerate() {
            assert_eq!(screen.menu_items[i].title, *expected_title);
            assert!(
                screen.menu_items[i].enabled,
                "Menu item '{}' should be enabled",
                expected_title
            );
        }
    }

    #[test]
    fn test_navigation_up() {
        let mut screen = create_test_screen();

        // Start at index 0, up should wrap to last item
        assert_eq!(screen.selected_item, 0);

        let key = KeyEvent::new(KeyCode::Up, KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_none()); // Should not change screen
        assert_eq!(screen.selected_item, screen.menu_items.len() - 1); // Wrapped to last
    }

    #[test]
    fn test_navigation_down() {
        let mut screen = create_test_screen();

        // Start at index 0, down should go to index 1
        assert_eq!(screen.selected_item, 0);

        let key = KeyEvent::new(KeyCode::Down, KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_none()); // Should not change screen
        assert_eq!(screen.selected_item, 1);
    }

    #[test]
    fn test_navigation_down_wrap() {
        let mut screen = create_test_screen();

        // Move to last item
        screen.selected_item = screen.menu_items.len() - 1;

        let key = KeyEvent::new(KeyCode::Down, KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_none()); // Should not change screen
        assert_eq!(screen.selected_item, 0); // Should wrap to first
    }

    #[test]
    fn test_enter_key_navigation() {
        let mut screen = create_test_screen();

        // Test selecting first item (File Selection)
        screen.selected_item = 0;
        let key = KeyEvent::new(KeyCode::Enter, KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_some());
        assert_eq!(result.unwrap(), ScreenEnum::FileSelection);
    }

    #[test]
    fn test_number_key_navigation() {
        let mut screen = create_test_screen();

        // Test pressing '2' should select Configuration and navigate there
        let key = KeyEvent::new(KeyCode::Char('2'), KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_some());
        assert_eq!(result.unwrap(), ScreenEnum::Configuration);
        assert_eq!(screen.selected_item, 1); // Should update selection too
    }

    #[test]
    fn test_number_key_out_of_range() {
        let mut screen = create_test_screen();

        // Test pressing '9' (out of range)
        let key = KeyEvent::new(KeyCode::Char('9'), KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_none()); // Should not navigate
        assert_eq!(screen.selected_item, 0); // Should not change selection
    }

    #[test]
    fn test_help_key() {
        let mut screen = create_test_screen();

        let key = KeyEvent::new(KeyCode::Char('h'), KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_some());
        assert_eq!(result.unwrap(), ScreenEnum::Help);
    }

    #[test]
    fn test_f1_help_key() {
        let mut screen = create_test_screen();

        let key = KeyEvent::new(KeyCode::F(1), KeyModifiers::empty());
        let result = screen.handle_key(key).unwrap();

        assert!(result.is_some());
        assert_eq!(result.unwrap(), ScreenEnum::Help);
    }

    #[test]
    fn test_all_menu_items_navigable() {
        let mut screen = create_test_screen();

        let expected_screens = vec![
            ScreenEnum::FileSelection,
            ScreenEnum::Configuration,
            ScreenEnum::Analysis,
            ScreenEnum::Database,
            ScreenEnum::Results,
            ScreenEnum::Help,
        ];

        for (i, expected_screen) in expected_screens.iter().enumerate() {
            screen.selected_item = i;
            let key = KeyEvent::new(KeyCode::Enter, KeyModifiers::empty());
            let result = screen.handle_key(key).unwrap();

            assert!(result.is_some(), "Menu item {} should be navigable", i);
            assert_eq!(
                result.unwrap(),
                *expected_screen,
                "Menu item {} should navigate to correct screen",
                i
            );
        }
    }

    #[test]
    fn test_unhandled_keys() {
        let mut screen = create_test_screen();

        // Test keys that should not do anything
        let unhandled_keys = vec![
            KeyCode::Tab,
            KeyCode::Backspace,
            KeyCode::Delete,
            KeyCode::Left,
            KeyCode::Right,
            KeyCode::Char('z'),
            KeyCode::F(5),
        ];

        for key_code in unhandled_keys {
            let key = KeyEvent::new(key_code, KeyModifiers::empty());
            let result = screen.handle_key(key).unwrap();

            assert!(
                result.is_none(),
                "Key {:?} should not trigger navigation",
                key_code
            );
        }
    }

    #[test]
    fn test_menu_item_descriptions() {
        let screen = create_test_screen();

        // Verify descriptions exist and are not empty
        for item in &screen.menu_items {
            assert!(
                !item.description.is_empty(),
                "Menu item '{}' should have a description",
                item.title
            );
            assert!(
                item.description.len() > 10,
                "Menu item '{}' should have a meaningful description",
                item.title
            );
        }
    }

    #[test]
    fn test_screen_enum_matches() {
        let screen = create_test_screen();

        // Verify that menu items map to correct screen enums
        let expected_mappings = vec![
            (0, ScreenEnum::FileSelection),
            (1, ScreenEnum::Configuration),
            (2, ScreenEnum::Analysis),
            (3, ScreenEnum::Database),
            (4, ScreenEnum::Results),
            (5, ScreenEnum::Help),
        ];

        for (index, expected_screen) in expected_mappings {
            assert_eq!(screen.menu_items[index].screen, expected_screen);
        }
    }

    #[test]
    fn test_case_insensitive_help() {
        let mut screen = create_test_screen();

        // Test both lowercase and uppercase 'h' for help
        let key_lower = KeyEvent::new(KeyCode::Char('h'), KeyModifiers::empty());
        let result_lower = screen.handle_key(key_lower).unwrap();

        let key_upper = KeyEvent::new(KeyCode::Char('H'), KeyModifiers::empty());
        let result_upper = screen.handle_key(key_upper).unwrap();

        assert_eq!(result_lower, result_upper);
        assert_eq!(result_lower, Some(ScreenEnum::Help));
    }
}
