Task: Refactor for clarity and maintainability with zero behavior change.

Rules:
- No public API changes, no semantic changes.
- Reduce duplication, improve names where ambiguity exists.
- Add comments where lifetimes/unsafe/Send+Sync bounds are subtle.

Acceptance:
- All tests pass unchanged.
- Clippy lint level: pedantic clean or justified allows.
