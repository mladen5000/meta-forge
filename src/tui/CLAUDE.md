# TUI Module Guidelines

## Purpose
Terminal user interface for interactive metagenomic analysis.

## Key Components
- `app.rs` - Main application state and event handling
- `components/` - Reusable UI components (navigation, status bar)
- `screens/` - Individual application screens
- `widgets/` - Custom UI widgets (progress bars, etc.)

## Development Rules
- Handle terminal resize events gracefully
- Use non-blocking I/O for responsiveness
- Provide keyboard shortcuts for all functions
- Test on different terminal emulators

## UI Design Principles
- Clear visual hierarchy and navigation
- Consistent color scheme and styling
- Accessible design with screen reader support
- Responsive layout for different terminal sizes

## State Management
- Use immutable state updates where possible
- Handle concurrent updates from background tasks
- Persist user preferences between sessions
- Validate user input before processing

## Performance Guidelines
- Minimize screen redraws for smooth interaction
- Use efficient data structures for large result sets
- Stream large outputs rather than loading all at once
- Profile rendering performance with large datasets

## Testing Approach
- Unit tests for UI components and state logic
- Integration tests for complete user workflows
- Test keyboard and mouse interaction handling
- Validate accessibility features

## Error Handling
- Display user-friendly error messages
- Provide recovery options for common errors
- Log detailed errors for debugging
- Handle terminal-related errors gracefully