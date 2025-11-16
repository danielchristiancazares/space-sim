# Repository Guidelines

## Project Structure & Module Organization
- `src/` contains the Bevy game logic; `src/atmosphere/` is the fluid core (grid, steps, diagnostics), and `src/player`, `src/tilemap`, etc. provide gameplay context.
- `tests/` hosts integration and regression suites that exercise multi-system flows; unit tests usually live beside the code they cover.
- `assets/` holds GLTF scenes, textures, and tuning CSVs used at runtime.
- `docs/` and top-level research logs capture design decisions—consult them before altering physics constants or solvers to avoid repeating past experiments.

## Build, Test, and Development Commands
- `cargo check` – fast structural verification; run before iterating on features.
- `cargo fmt && cargo clippy --all-targets -- -D warnings` – enforces style and lints; CI treats warnings as errors.
- `cargo test` – executes unit and integration tests under both `src/**` and `tests/`.
- `cargo run --features bevy_winit` (default) or `cargo run --release` – launches the simulation client; use `--release` when profiling atmosphere performance.

## Coding Style & Naming Conventions
- Rustfmt defaults (4-space indent, trailing commas) are mandatory; never hand-format files.
- Modules and functions use `snake_case`; types and systems use `UpperCamelCase`; constants reside in `constants.rs` and use `SCREAMING_SNAKE_CASE`.
- Keep atmospheric parameters centralized in `src/atmosphere/constants.rs`; cite doc links (e.g., `docs/bug_mass_conservation_leak.md`) when changing them.
- Prefer small, single-purpose systems wired through Bevy sets; favor explicit buffer resources over implicit global state.

## Testing Guidelines
- Use `#[cfg(test)]` blocks next to the logic for solver micro-tests; integration scenarios belong in `tests/`.
- Name tests after the behavior under scrutiny (`test_pressure_projection_reduces_divergence`) and assert on physical invariants (mass, momentum, positivity).
- When adding a new numerical method, include a regression test that reproduces the motivating scenario plus a diagnostic snapshot in `docs/`.
- Run `cargo test -- --nocapture` when investigating divergence logs so failures surface detailed traces.

## Commit & Pull Request Guidelines
- Commits follow imperative, sentence-case summaries (“Fix mass tracker drift”) with optional extended descriptions explaining physics rationale.
- Each PR should: describe the scenario reproduced, list verification commands (`cargo test`, targeted Bevy runs), reference relevant docs/issues, and attach screenshots or CSV snippets when visual/diagnostic output changes.
- Highlight expected impacts on stability, performance, and tuning knobs so reviewers can compare with historical log entries under `docs/`.
