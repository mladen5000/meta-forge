# Text User Interface (TUI) Architecture Design
## Metagenomic Analysis Pipeline

### Executive Summary

This document outlines the architecture for a comprehensive Text User Interface (TUI) that will replace the current command-line interface for the metagenomic analysis pipeline. The TUI will provide an intuitive, real-time, and interactive experience for both expert bioinformaticians and newcomers to metagenomic analysis.

---

## 1. Overall TUI Architecture & Component Hierarchy

### 1.1 Core Architecture Principles

```
┌─────────────────────────────────────────────────────────────────┐
│                        TUI Application                          │
├─────────────────────────────────────────────────────────────────┤
│  App State Manager (Global State)                              │
│  ├── Pipeline State                                            │
│  ├── UI State                                                  │
│  ├── Configuration State                                       │
│  └── Progress State                                            │
├─────────────────────────────────────────────────────────────────┤
│  Event System (Centralized Event Handling)                     │
│  ├── Keyboard Events                                           │
│  ├── Pipeline Events                                           │
│  ├── Timer Events                                              │
│  └── Database Events                                           │
├─────────────────────────────────────────────────────────────────┤
│  Screen Manager (Navigation & Layout)                          │
│  ├── Dashboard Screen                                          │
│  ├── Analysis Screen                                           │
│  ├── Database Screen                                           │
│  ├── Configuration Screen                                      │
│  └── Results Screen                                            │
├─────────────────────────────────────────────────────────────────┤
│  Widget Library (Reusable Components)                          │
│  ├── File Browser                                              │
│  ├── Progress Visualizations                                   │
│  ├── Interactive Forms                                         │
│  ├── Data Tables                                               │
│  └── Status Indicators                                         │
├─────────────────────────────────────────────────────────────────┤
│  Pipeline Integration Layer                                     │
│  ├── Command Adapter                                           │
│  ├── Progress Monitor                                          │
│  ├── Result Parser                                             │
│  └── Configuration Bridge                                      │
└─────────────────────────────────────────────────────────────────┘
```

### 1.2 Technology Stack

- **Core TUI Framework**: [ratatui](https://github.com/ratatui-org/ratatui) v0.24+
- **Terminal Backend**: [crossterm](https://github.com/crossterm-rs/crossterm) v0.27+
- **State Management**: Custom reactive state system
- **Event Handling**: Tokio-based async event loop
- **File Operations**: Built-in Rust std + custom file browser
- **Progress Tracking**: Real-time pipeline integration

---

## 2. Screen Layouts & Navigation Flow

### 2.1 Main Dashboard Screen

```
┌─ MetaForge Dashboard ──────────────────────────────────────────────────────────┐
│ [F1: Help] [F2: Config] [F3: Database] [F4: Analysis] [F5: Results] [Q: Quit]  │
├────────────────────────────────────────────────────────────────────────────────┤
│ ┌─ System Status ────────────────┐ ┌─ Recent Analysis ─────────────────────────┐ │
│ │ ● Database: Connected          │ │ Sample_001    | Complete  | 2h 15m ago   │ │
│ │ ● Memory: 2.3GB / 16GB        │ │ Sample_002    | Running   | Progress 45% │ │
│ │ ● CPU: 45% (8 cores)          │ │ Sample_003    | Failed    | 1d ago       │ │
│ │ ● Temp Space: 15GB available  │ │ Sample_004    | Queued    | -            │ │
│ └────────────────────────────────┘ └───────────────────────────────────────────┘ │
│ ┌─ Quick Actions ─────────────────────────────────────────────────────────────┐ │
│ │ [N] New Analysis   [D] Database Manager   [C] Configuration   [R] Results  │ │
│ └─────────────────────────────────────────────────────────────────────────────┘ │
│ ┌─ Pipeline Overview ─────────────────────────────────────────────────────────┐ │
│ │ Active Processes: 2                                                         │ │
│ │                                                                             │ │
│ │ 1. Sample_002_analysis                                                      │ │
│ │    ├─ Read Processing    ████████████████████░░░░░  75% (2.1M/2.8M reads) │ │
│ │    ├─ Assembly          ██████████░░░░░░░░░░░░░░░░  40% (Building graph)  │ │
│ │    └─ Feature Extract   ░░░░░░░░░░░░░░░░░░░░░░░░░░   0% (Pending)         │ │
│ │                                                                             │ │
│ │ 2. Database_indexing                                                        │ │
│ │    └─ K-mer Index       ██████████████████████████ 100% (Complete)        │ │
│ └─────────────────────────────────────────────────────────────────────────────┘ │
│ ┌─ Logs ──────────────────────────────────────────────────────────────────────┐ │
│ │ [INFO ] Starting analysis for Sample_002                                    │ │
│ │ [DEBUG] Using 8 threads for parallel processing                            │ │
│ │ [WARN ] Low memory available, reducing cache size                          │ │
│ │ [INFO ] Assembly graph construction completed (N50: 15,432 bp)             │ │
│ └─────────────────────────────────────────────────────────────────────────────┘ │
└────────────────────────────────────────────────────────────────────────────────┘
```

### 2.2 Analysis Configuration Screen

```
┌─ New Analysis Configuration ─────────────────────────────────────────────────────┐
│ [Tab: Next Field] [Enter: Confirm] [Esc: Cancel]                                │
├──────────────────────────────────────────────────────────────────────────────────┤
│ ┌─ Input Files ─────────────────────┐ ┌─ Analysis Parameters ──────────────────┐ │
│ │ Sample Name: [Sample_005_____] *  │ │ Analysis Mode:                          │ │
│ │                                   │ │ ○ Fast (Basic features, ~30 min)       │ │
│ │ Input Type:                       │ │ ● Standard (Full pipeline, ~2h)        │ │
│ │ ● Paired-end FASTQ               │ │ ○ Accurate (Advanced ML, ~6h)          │ │
│ │ ○ Single-end FASTQ               │ │ ○ Custom                                │ │
│ │ ○ Assembled contigs               │ │                                         │ │
│ │                                   │ │ Threads: [8_] (1-16 available)         │ │
│ │ Files Selected:                   │ │ Memory Limit: [8_] GB (max 16 GB)      │ │
│ │ □ /data/R1.fastq (2.1 GB)        │ │                                         │ │
│ │ □ /data/R2.fastq (2.0 GB)        │ │ K-mer Range: [21] to [41]               │ │
│ │                                   │ │ Min Coverage: [2_]                      │ │
│ │ [Browse Files...]                 │ │                                         │ │
│ └───────────────────────────────────┘ └─────────────────────────────────────────┘ │
│ ┌─ Output Configuration ────────────────────────────────────────────────────────┐ │
│ │ Output Directory: [/results/Sample_005_________________] [Browse...]          │ │
│ │                                                                              │ │
│ │ Output Formats:                                                              │ │
│ │ ☑ JSON (Machine readable)    ☑ CSV (Spreadsheet compatible)                │ │
│ │ ☑ FASTA (Sequences)          □ Binary (Fast loading)                       │ │
│ │                                                                              │ │
│ │ Include in Results:                                                          │ │
│ │ ☑ Assembly contigs          ☑ Taxonomic classification                      │ │
│ │ ☑ Feature vectors          ☑ Abundance profiles                            │ │
│ │ ☑ Quality metrics          □ Debug information                             │ │
│ └──────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ [Start Analysis] [Save as Template] [Load Template] [Preview Command]          │ └──────────────────────────────────────────────────────────────────────────────────┘
```

### 2.3 File Browser Widget

```
┌─ Select Input Files ─────────────────────────────────────────────────────────────┐
│ Path: /data/metagenomics/samples/                         [↑] [Home] [Refresh]  │
├──────────────────────────────────────────────────────────────────────────────────┤
│ Directory:                                                                       │
│ drwxr-xr-x  4096  ../                                                           │
│ drwxr-xr-x  4096  sample_001/                                                   │ 
│ drwxr-xr-x  4096  sample_002/                                                   │
│ drwxr-xr-x  4096  sample_003/                                                   │
│ -rw-r--r--  2.1G  ● SRR390728_1.fastq    [FASTQ, Paired R1]                   │
│ -rw-r--r--  2.0G  ● SRR390728_2.fastq    [FASTQ, Paired R2]                   │
│ -rw-r--r--  156K  readme.txt                                                    │
│ -rw-r--r--  1.2M  metadata.csv                                                  │
│                                                                                  │
│ File Details:                                                                    │
│ ┌─ SRR390728_1.fastq ────────────────────────────────────────────────────────┐ │
│ │ Size: 2,147,483,648 bytes (2.1 GB)                                         │ │
│ │ Type: FASTQ sequence file                                                   │ │
│ │ Reads: ~8,500,000 (estimated)                                              │ │
│ │ Quality: Illumina 1.8+ encoding                                            │ │
│ │ Paired: Yes (R1, pair with SRR390728_2.fastq)                              │ │
│ └─────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ Selected Files: SRR390728_1.fastq, SRR390728_2.fastq                           │
│ [Select] [Select All] [Clear] [Cancel]                                         │
└──────────────────────────────────────────────────────────────────────────────────┘
```

### 2.4 Real-time Progress Screen

```
┌─ Analysis Progress: Sample_005 ──────────────────────────────────────────────────┐
│ Started: 2024-01-15 14:30:22 | Elapsed: 01:23:45 | Est. Remaining: 00:45:12     │
├──────────────────────────────────────────────────────────────────────────────────┤
│ ┌─ Pipeline Stages ────────────────────────────────────────────────────────────┐ │
│ │ ✓ Input Validation        00:00:05   [Complete]                              │ │
│ │ ✓ Quality Control        00:12:30   [Complete - 95% reads passed]           │ │
│ │ ▶ Read Processing        00:45:22   [Running - 8.5M/10.2M reads]            │ │
│ │   ████████████████████████████████████████████████████████████░░░░░ 83%     │ │
│ │ ⏸ Assembly Construction  --:--:--   [Pending]                               │ │
│ │ ⏸ Feature Extraction     --:--:--   [Pending]                               │ │
│ │ ⏸ Taxonomic Classification --:--:--  [Pending]                              │ │
│ │ ⏸ Report Generation      --:--:--   [Pending]                               │ │
│ └──────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ ┌─ Current Stage Details ──────────────────────────────────────────────────────┐ │
│ │ Stage: Read Processing & Quality Filtering                                   │ │
│ │                                                                              │ │
│ │ Throughput: 185,432 reads/sec                                               │ │
│ │ Memory Usage: 5.2 GB / 8.0 GB allocated                                     │ │
│ │ CPU Usage: 94% (8/8 cores active)                                           │ │
│ │ I/O Rate: Read 145 MB/s, Write 23 MB/s                                      │ │
│ │                                                                              │ │
│ │ Quality Metrics:                                                             │ │
│ │ • Reads processed: 8,543,221 / 10,200,000                                   │ │
│ │ • Reads passed QC: 8,121,459 (95.1%)                                        │ │
│ │ • Reads discarded: 421,762 (4.9%)                                           │ │
│ │ • Mean read length: 142 bp                                                  │ │
│ │ • Mean quality score: 35.2                                                  │ │
│ └──────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ ┌─ System Resources ───────────────────────────────────────────────────────────┐ │
│ │ CPU:  ████████████████████████████████████████████████████░░  94% (8 cores)  │ │
│ │ MEM:  ████████████████████████████████████░░░░░░░░░░░░░░░░░░ 65% (10.4/16 GB)│ │
│ │ DISK: ██████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 12% (89/750 GB) │ │
│ │ TEMP: ████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 23% (17/75 GB) │ │
│ └──────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ [Pause] [Stop] [View Logs] [Advanced Details]                                  │
└──────────────────────────────────────────────────────────────────────────────────┘
```

### 2.5 Database Management Screen

```
┌─ Database Manager ───────────────────────────────────────────────────────────────┐
│ Connected: /data/metagenomics.db (157 MB) | Status: ● Active                     │
├──────────────────────────────────────────────────────────────────────────────────┤
│ ┌─ Database Operations ─────────────┐ ┌─ Database Contents ──────────────────────┐ │
│ │                                   │ │ Tables:                                  │ │
│ │ [Initialize New Database]         │ │ • samples (1,247 records)               │ │
│ │                                   │ │ • sequences (156,432 records)           │ │
│ │ [Import Taxonomy Data]            │ │ • taxonomy (45,123 records)             │ │
│ │                                   │ │ • k_mers (2,847,293 records)            │ │
│ │ [Import Reference Sequences]      │ │ • analysis_results (1,247 records)      │ │
│ │                                   │ │                                          │ │
│ │ [Build K-mer Index]              │ │ Indices:                                 │ │
│ │                                   │ │ ✓ Primary keys                          │ │
│ │ [Vacuum & Optimize]               │ │ ✓ K-mer hash index (97% coverage)       │ │
│ │                                   │ │ ✓ Taxonomy tree index                   │ │
│ │ [Export Data]                     │ │ ⚠ Sequence similarity index (missing)   │ │
│ │                                   │ │                                          │ │
│ │ [Database Statistics]             │ │ Database Health: ● Good                  │ │
│ │                                   │ │ Last vacuum: 2 days ago                 │ │
│ │ [Connection Settings]             │ │ Fragmentation: 12%                       │ │
│ └───────────────────────────────────┘ └──────────────────────────────────────────┘ │
│                                                                                  │ │
│ ┌─ Recent Samples ─────────────────────────────────────────────────────────────┐ │
│ │ ID    | Sample Name    | Date       | Status    | Contigs | Features         │ │
│ │ 1247  | Sample_005     | 2024-01-15 | Complete  | 1,234   | ✓               │ │
│ │ 1246  | Sample_004     | 2024-01-14 | Complete  | 2,156   | ✓               │ │
│ │ 1245  | Sample_003     | 2024-01-14 | Failed    | 0       | ✗               │ │
│ │ 1244  | Sample_002     | 2024-01-13 | Complete  | 3,421   | ✓               │ │
│ │ 1243  | Sample_001     | 2024-01-13 | Complete  | 1,876   | ✓               │ │
│ └──────────────────────────────────────────────────────────────────────────────┘ │
│                                                                                  │
│ [Back to Dashboard] [Query Builder] [Export Results]                            │
└──────────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Widget Selection & Placement Strategy

### 3.1 Widget Hierarchy

```rust
// Core widget types and their purposes
pub enum WidgetType {
    // Layout widgets
    Container,      // Main layout containers
    Panel,          // Bordered sections
    Tabs,           // Multi-page content
    
    // Input widgets
    TextInput,      // Text entry fields
    NumberInput,    // Numeric entry with validation
    FileSelect,     // File/directory selection
    Checkbox,       // Boolean options
    RadioGroup,     // Mutually exclusive options
    
    // Display widgets
    ProgressBar,    // Linear progress indicators
    ProgressGauge,  // Circular/gauge progress
    Table,          // Tabular data display
    Tree,           // Hierarchical data
    Chart,          // Data visualization
    LogViewer,      // Scrollable log display
    
    // Interactive widgets
    Menu,           // Navigation menus
    Button,         // Action triggers
    StatusBar,      // System status display
}
```

### 3.2 Widget Placement Strategy

#### 3.2.1 Layout Grid System
- **12-column responsive grid** for consistent alignment
- **Minimum widget sizes** to ensure readability
- **Flexible containers** that adapt to terminal size
- **Priority-based overflow handling** for small terminals

#### 3.2.2 Widget Priority Classes
1. **Critical**: Always visible (status bar, navigation)
2. **Primary**: Core functionality (main panels)
3. **Secondary**: Supporting information (details, logs)
4. **Optional**: Nice-to-have features (charts, advanced options)

#### 3.2.3 Responsive Behavior
```rust
pub struct ResponsiveLayout {
    // Terminal size breakpoints
    small: (u16, u16),   // < 80x24
    medium: (u16, u16),  // 80x24 to 120x40
    large: (u16, u16),   // > 120x40
    
    // Widget visibility rules per breakpoint
    widget_rules: HashMap<WidgetType, VisibilityRule>,
}
```

---

## 4. Integration Points with Existing Pipeline Code

### 4.1 Command Adapter Pattern

```rust
pub struct CommandAdapter {
    /// Converts TUI interactions to CLI commands
    pub fn adapt_analysis_request(&self, ui_request: AnalysisRequest) -> Vec<String> {
        // Convert UI parameters to CLI arguments
        // Handle file paths, options, and configurations
        // Return command strings for pipeline execution
    }
    
    pub fn adapt_database_operation(&self, ui_request: DatabaseRequest) -> DatabaseOps {
        // Convert UI database operations to internal structures
    }
}
```

### 4.2 Progress Monitor Bridge

```rust
pub struct ProgressMonitor {
    /// Integrate with existing progress_display.rs
    pipeline_progress: Arc<RwLock<PipelineProgress>>,
    
    pub fn subscribe_to_pipeline(&mut self) -> mpsc::Receiver<ProgressUpdate> {
        // Connect to existing progress tracking system
    }
    
    pub fn convert_progress(&self, update: ProgressUpdate) -> TUIProgressUpdate {
        // Convert pipeline progress to TUI-compatible format
    }
}
```

### 4.3 Configuration Bridge

```rust
pub struct ConfigurationBridge {
    /// Bridge TUI settings to PipelineConfiguration
    pub fn sync_ui_to_pipeline(&self, ui_config: TUIConfig) -> PipelineConfiguration {
        // Convert TUI configuration to pipeline configuration
    }
    
    pub fn sync_pipeline_to_ui(&self, pipeline_config: PipelineConfiguration) -> TUIConfig {
        // Convert pipeline configuration to TUI format
    }
}
```

### 4.4 Result Parser Integration

```rust
pub struct ResultParser {
    /// Parse pipeline outputs for TUI display
    pub fn parse_analysis_results(&self, results: AnalysisResults) -> TUIResultView {
        // Convert complex pipeline results to TUI-friendly format
        // Handle visualizations, tables, and summaries
    }
}
```

---

## 5. State Management Approach

### 5.1 Global State Architecture

```rust
pub struct AppState {
    // Core application state
    current_screen: Screen,
    previous_screen: Vec<Screen>, // Navigation history
    
    // Pipeline integration state
    pipeline_state: PipelineState,
    progress_state: ProgressState,
    
    // UI-specific state
    ui_state: UIState,
    modal_stack: Vec<Modal>,
    
    // Configuration state
    config: TUIConfiguration,
    pipeline_config: PipelineConfiguration,
}

#[derive(Clone)]
pub enum PipelineState {
    Idle,
    Running {
        sample_name: String,
        stage: PipelineStage,
        progress: f64,
        started_at: Instant,
    },
    Completed {
        results: AnalysisResults,
        duration: Duration,
    },
    Failed {
        error: String,
        stage: PipelineStage,
    },
}
```

### 5.2 Reactive State Updates

```rust
pub trait StateSubscriber {
    fn on_state_change(&mut self, change: StateChange);
}

pub struct StateManager {
    state: Arc<RwLock<AppState>>,
    subscribers: Vec<Box<dyn StateSubscriber + Send + Sync>>,
    
    pub fn subscribe(&mut self, subscriber: Box<dyn StateSubscriber + Send + Sync>) {
        self.subscribers.push(subscriber);
    }
    
    pub fn update_state<F>(&self, updater: F) 
    where F: FnOnce(&mut AppState) {
        let mut state = self.state.write().unwrap();
        updater(&mut state);
        // Notify all subscribers of changes
        self.notify_subscribers(&state);
    }
}
```

---

## 6. Event Handling System

### 6.1 Event Types

```rust
pub enum TUIEvent {
    // User input events
    KeyPress(KeyEvent),
    Mouse(MouseEvent),
    
    // System events
    Terminal(TerminalEvent),
    Timer(TimerEvent),
    
    // Pipeline events
    Pipeline(PipelineEvent),
    Progress(ProgressEvent),
    
    // Application events
    Navigation(NavigationEvent),
    Modal(ModalEvent),
}

pub enum PipelineEvent {
    Started { sample_name: String },
    StageChanged { stage: PipelineStage },
    Progress { stage: PipelineStage, progress: f64 },
    Completed { results: AnalysisResults },
    Failed { error: String, stage: PipelineStage },
    Log { level: LogLevel, message: String },
}
```

### 6.2 Event Loop Architecture

```rust
pub struct EventLoop {
    event_tx: mpsc::Sender<TUIEvent>,
    event_rx: mpsc::Receiver<TUIEvent>,
    
    // Event handlers
    handlers: HashMap<EventType, Box<dyn EventHandler>>,
    
    pub async fn run(&mut self) -> Result<()> {
        loop {
            tokio::select! {
                // Handle keyboard/terminal events
                event = self.terminal_events() => {
                    self.handle_terminal_event(event?).await?;
                }
                
                // Handle pipeline events
                event = self.pipeline_events() => {
                    self.handle_pipeline_event(event?).await?;
                }
                
                // Handle timer events (for UI updates)
                _ = self.ui_timer.tick() => {
                    self.handle_ui_tick().await?;
                }
            }
        }
    }
}
```

### 6.3 Keyboard Shortcut System

```rust
pub struct KeybindManager {
    bindings: HashMap<KeyCombo, Action>,
    context_bindings: HashMap<Screen, HashMap<KeyCombo, Action>>,
    
    pub fn register_global_binding(&mut self, key: KeyCombo, action: Action) {
        self.bindings.insert(key, action);
    }
    
    pub fn register_context_binding(&mut self, screen: Screen, key: KeyCombo, action: Action) {
        self.context_bindings.entry(screen).or_insert_with(HashMap::new)
            .insert(key, action);
    }
}

// Global shortcuts
// F1: Help, F2: Config, F3: Database, F4: Analysis, F5: Results
// Ctrl+Q: Quit, Ctrl+C: Cancel current operation
// Tab/Shift+Tab: Navigate between fields
// Enter: Confirm/Execute, Esc: Cancel/Back
```

---

## 7. Architecture Decision Records (ADRs)

### 7.1 ADR-001: TUI Framework Selection

**Status**: Accepted  
**Decision**: Use ratatui + crossterm for TUI implementation  
**Rationale**: 
- ratatui provides excellent widget ecosystem and layout management
- crossterm offers reliable cross-platform terminal handling
- Strong Rust ecosystem integration
- Active maintenance and community support

**Alternatives Considered**:
- cursive: Less flexible layout system
- tui-rs: Deprecated in favor of ratatui
- Custom terminal handling: Too much low-level complexity

### 7.2 ADR-002: State Management Pattern

**Status**: Accepted  
**Decision**: Centralized state with reactive updates  
**Rationale**:
- Single source of truth for application state
- Predictable state transitions
- Easy integration with existing pipeline code
- Supports real-time updates from background processes

**Trade-offs**:
- Slightly more complex than local component state
- Requires careful state synchronization
- Benefits outweigh complexity for this use case

### 7.3 ADR-003: Pipeline Integration Strategy

**Status**: Accepted  
**Decision**: Adapter pattern for CLI command translation  
**Rationale**:
- Minimal changes to existing pipeline code
- Clear separation of concerns
- Easy testing and maintenance
- Preserves existing functionality

### 7.4 ADR-004: Progress Visualization Approach

**Status**: Accepted  
**Decision**: Multi-level progress display with real-time updates  
**Rationale**:
- Better user experience for long-running operations
- Provides meaningful feedback at multiple granularities
- Integrates well with existing progress tracking
- Supports both determinate and indeterminate progress

---

## 8. Implementation Phases

### Phase 1: Core Infrastructure (Weeks 1-2)
- Basic TUI application structure
- Event system implementation
- State management foundation
- Basic navigation between screens

### Phase 2: Essential Screens (Weeks 3-4)
- Dashboard screen implementation
- Analysis configuration screen
- File browser widget
- Basic progress display

### Phase 3: Pipeline Integration (Weeks 5-6)
- Command adapter implementation
- Progress monitor integration
- Configuration bridge
- Basic analysis execution

### Phase 4: Advanced Features (Weeks 7-8)
- Database management screen
- Results visualization
- Advanced progress displays
- Error handling and recovery

### Phase 5: Polish & Testing (Weeks 9-10)
- Comprehensive testing
- Performance optimization
- Documentation
- User experience refinements

---

## 9. Quality Attributes & Non-Functional Requirements

### 9.1 Performance Requirements
- **Startup Time**: < 1 second on modern hardware
- **UI Responsiveness**: < 100ms response to user input
- **Memory Usage**: < 50MB for TUI interface (excluding pipeline)
- **Terminal Compatibility**: Support terminals ≥ 80x24 characters

### 9.2 Usability Requirements
- **Learning Curve**: New users productive within 15 minutes
- **Help System**: Context-sensitive help available on F1
- **Error Recovery**: Graceful handling of pipeline failures
- **Accessibility**: Screen reader compatible (where possible)

### 9.3 Reliability Requirements
- **Pipeline Integration**: No data loss during TUI crashes
- **State Persistence**: Resume operations after restart
- **Error Handling**: Clear error messages with suggested actions
- **Resource Monitoring**: Prevent system resource exhaustion

---

## 10. Future Extensions

### 10.1 Advanced Visualizations
- Real-time assembly graph visualization
- Interactive taxonomic trees
- Quality score distributions
- Coverage plots

### 10.2 Workflow Management
- Save/load analysis templates
- Batch processing interface
- Pipeline scheduling
- Results comparison tools

### 10.3 Remote Operation
- Remote pipeline execution
- Distributed computing support
- Cloud storage integration
- Collaborative analysis features

---

## Conclusion

This TUI architecture provides a comprehensive, user-friendly interface that transforms the complex metagenomic analysis pipeline into an accessible, interactive experience. The design prioritizes usability while maintaining the powerful functionality of the existing CLI tools, providing both novice and expert users with an intuitive path to genomic insights.

The modular architecture ensures maintainability and extensibility, while the integration strategy preserves existing investments in the pipeline codebase. The phased implementation approach allows for iterative development and early user feedback.