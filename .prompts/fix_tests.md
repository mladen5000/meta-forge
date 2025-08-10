Task: Make the Rust test suite pass with minimal changes.

Context:
- Failing tests and compiler/clippy output are included below.
- Do not change public APIs or feature flags unless strictly required by the failure.

Acceptance:
- `make check` succeeds (fmt, clippy, check, tests).
- Diff is surgical.
- If a real bug is fixed, add/adjust tests demonstrating the bug and the fix.
