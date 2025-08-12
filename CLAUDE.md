# CLAUDE.md – Project Guidelines for Claude Code

Claude, follow these rules when working on this project.

---

## 🎯 Primary Goals
- Maintain a **clippy-clean**, idiomatic Rust codebase.
- Keep diffs **surgical** — only touch what’s needed.
- Preserve public APIs unless explicitly told otherwise.
- Optimize for readability, correctness, and maintainability.

---

## ⚙️ Development Workflow

1. **Plan before coding**
   - Read the relevant files and dependencies.
   - Outline a step-by-step plan in the scratchpad.
   - Confirm the plan before implementing.

2. **Test-Driven Development (TDD)**
   - Write failing tests first.
   - Implement the minimal code to make them pass.
   - Refactor with tests still passing.

3. **Small Commits**
   - Make focused changes in separate commits.
   - Each commit must pass `make check`.

4. **Code Review Mindset**
   - Add inline comments where lifetimes, `unsafe`, or concurrency behavior may be non-obvious.
   - Flag any potential performance bottlenecks.

---

## 🛠 Commands You May Run
- **Formatting & Linting**
  - `cargo fmt --all`
  - `cargo clippy --all-targets --all-features -- # -D warnings`
- **Compilation & Checks**
  <!-- - `cargo check --all-targets --all-features` -->
- **Testing**
  <!-- - `cargo test` (or `cargo nextest run` if available)
  - `cargo test -- --quiet --skip slow` (smoke tests) -->
- **Benchmarking**
  - `cargo bench` (only when asked)

---

## 📦 Code Style
- Edition: Rust 2021 or 2024 - choose 1 and stick with it.
- Prefer explicit over implicit lifetimes if clarity improves.
- Use `?` for error propagation where possible.
- Avoid panics in library code; use `Result` or custom errors.
- Keep functions short and focused.

---

## 🚫 Things You May Not Do
- Do not change public API without request.
- Do not run arbitrary network calls.
- Do not reformat unrelated files.
- Do not bypass safety checks in `unsafe` code without a clear justification in comments.

---



## 🚨 CRITICAL: CONCURRENT EXECUTION & FILE MANAGEMENT

**ABSOLUTE RULES**:
1. ALL operations MUST be concurrent/parallel in a single message
2. **NEVER save working files, text/mds and tests to the root folder**
3. ALWAYS organize files in appropriate subdirectories

---

### ⚡ GOLDEN RULE: "1 MESSAGE = ALL RELATED OPERATIONS"

**MANDATORY PATTERNS:**
- **TodoWrite**: ALWAYS batch ALL todos in ONE call (5-10+ todos minimum)
- **Task tool**: ALWAYS spawn ALL agents in ONE message with full instructions
- **File operations**: ALWAYS batch ALL reads/writes/edits in ONE message
- **Bash commands**: ALWAYS batch ALL terminal operations in ONE message
- **Memory operations**: ALWAYS batch ALL memory store/retrieve in ONE message

---

### 📁 File Organization Rules

**NEVER save to root folder. Use these directories:**
- `/src` - Source code files
- `/tests` - Test files
- `/docs` - Documentation and markdown files
- `/config` - Configuration files
- `/scripts` - Utility scripts
- `/examples` - Example code
- `/benches` - Benchmarks


---

## Code Style & Best Practices

- **Modular Design**: Files under 500 lines
- **Environment Safety**: Never hardcode secrets
- **Test-First**: Write tests before implementation
- **Clean Architecture**: Separate concerns
- **Documentation**: Keep updated

### Migration & Planning
`rust-bio-optimizer`, `agent-organizer`, `bio-testing-validation`

## 🎯 Claude Code vs MCP Tools

### Claude Code Handles ALL:
- File operations (Read, Write, Edit, MultiEdit, Glob, Grep)
- Code generation and programming
- Bash commands and system operations
- Implementation work
- Project navigation and analysis
- TodoWrite and task management
<!-- - Git operations -->
- Package management
- Testing and debugging

### MCP Tools ONLY:
- Coordination and planning
- Memory management
- Neural features
- Performance tracking
- Swarm orchestration
- GitHub integration

**KEY**: MCP coordinates, Claude Code executes.

## important-instruction-reminders
Do what has been asked; nothing more, nothing less.
NEVER create files unless they're absolutely necessary for achieving your goal.
ALWAYS prefer editing an existing file to creating a new one.
NEVER proactively create documentation files (*.md) or README files. Only create documentation files if explicitly requested by the User.
Never save working files, text/mds and tests to the root folder.

## Project Status

- Metagenomic assembly pipeline successfully implemented with TUI interface
- K-mer processing and graph construction optimized for performance
- Error correction and quality assessment working with ML enhancements
- Output generation and reporting complete with multiple formats
- SQLite database integration for taxonomy and k-mer storage
- Neural network models for repeat resolution and learned data structures

## Module Organization

Each `src/` subdirectory contains its own `CLAUDE.md` with specific guidelines:
- `src/assembly/` - Sequence assembly with adaptive k-mer strategies
- `src/core/` - Fundamental data structures and utilities
- `src/database/` - SQLite-based storage and querying
- `src/features/` - Feature extraction from sequences and graphs
- `src/ml/` - Machine learning models and neural networks
- `src/pipeline/` - High-level workflow orchestration
- `src/tui/` - Terminal user interface components
- `src/utils/` - Shared utilities and helper functions

## Project Memories
- Use the bioinformatics agents when possible
- Each module has specific performance and testing requirements
- Follow modular design principles for maintainability