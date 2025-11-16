# Space-Sim Atmosphere Prototype

This repo contains a Bevy-based research prototype that simulates a 2D cross-section of a space-station atmosphere. The focus is rapid iteration on grid-based fluid steps (advection, diffusion, compression heating) plus supporting systems such as life-support sources, respiration sinks, diagnostics, and CSV logging. Earlier pressure-flux and Poisson projection experiments were removed from the live loop, so the current build reflects only the components present under `src/atmosphere/simulation.rs`.

## Project Layout

| Path | Purpose |
| --- | --- |
| `src/main.rs`, `src/lib.rs` | Game entry point and Bevy app wiring |
| `src/atmosphere/` | Grid, constants, diagnostics, monitoring, CSV export, and the active steps (advection, diffusion, compression) |
| `src/tilemap.rs`, `src/player.rs`, `src/animation.rs` | World representation, player hooks, camera/UX logic |
| `tests/` | Integration/regression tests for atmosphere subsystems |
| `assets/` | Scenes, textures, and tuning data loaded at runtime |
| `docs/` | Research logs explaining removed pressure-flux/projection work and known limitations |

## Feature Overview

### Rendering Engine Systems
- **8-directional character movement** with normalized diagonal speed and per-axis collision
- **Smooth camera following** with configurable interpolation
- **Tile-based collision detection** using decoupled collision map
- **Character sprite rotation** for 8 compass directions
- **Real-time visualization** (press V): Pressure, temperature, O₂, N₂, CO₂, breathability, and velocity fields

## Current Simulation Pipeline

Each frame is subdivided by CFL constraints (`simulate_atmosphere`):

1. **Flux-Conservative Advection** (`steps/advection.rs`) – mass-conserving density transport plus MacCormack velocity/temperature advection.
2. **Explicit Diffusion** (`steps/diffusion.rs`) – viscosity, thermal diffusion, and scalar mixing, respecting tile collisions.
3. **Compression Heating** (`steps/compression.rs`) – applies ∂T/∂t = -(γ-1)T∇·u using the current divergence field.
4. **Pressure Update** (`simulation::update_grid_pressures`) – recomputes ideal-gas pressures for diagnostics and source terms.

Pressure-driven flux and Poisson pressure projection are *not* wired in; any mention of them in historical docs refers to deprecated code paths. Mass equalization presently occurs only through advection/diffusion and source terms, so expect lingering pressure gradients in longer runs.

## Build & Run

```bash
cargo check            # Fast type/build check
cargo test             # Run unit + integration tests (including diagnostics)
cargo fmt
cargo clippy --all-targets -- -D warnings
cargo run              # Launches the simulation (use --release when profiling)
```

The project targets Rust 1.75+ and Bevy 0.17. Assets assume X11 by default (`bevy_winit` feature is enabled in Cargo.toml).

## Testing & Diagnostics

- Unit tests live next to their modules; integration suites run from `tests/`.
- Diagnostics (`src/atmosphere/diagnostics.rs`) track mass, divergence, and conservation errors; CSV exports land under `.debug/diagnostics/`.
- When validating numerical changes, attach relevant CSV snippets or log excerpts to commits/PRs so reviewers can compare against prior runs.

### Atmospheric Simulation Pipeline
The atmosphere plugin orchestrates a four-step CFD loop each frame:

1. **Life support generation** (`sources.rs`): Adds O₂/N₂ at life support tiles based on `LIFE_SUPPORT_O2_RATE` (4.7 kg/s) and `LIFE_SUPPORT_N2_RATE` (15.6 kg/s)

2. **Forced convection** (`sources.rs`): Fans create 1.5 m/s tangential swirl within a 3-tile radius to distribute gases faster than diffusion alone

3. **Player respiration** (`sources.rs`): Consumes O₂ at 5.54×10⁻⁶ kg/s and produces CO₂ at 6.09×10⁻⁶ kg/s (human resting rate)

4. **Pressure update** (`simulation.rs`): Recomputes pressure from densities using the ideal gas law

5. **Fluid dynamics simulation** (`steps/`):
   - **Advection** (`steps/advection.rs`): MacCormack transport of density, momentum, and temperature
   - **Diffusion** (`steps/diffusion.rs`): Viscous and scalar diffusion
   - **Compression Heating** (`steps/compression.rs`): Temperature changes from velocity divergence

### Collision System
The `TileCollisionMap` (src/tilemap.rs:14-68) separates rendering from gameplay:
- Mirrors tilemap transform (origin, tile size) to convert world coordinates to tile indices
- Queried by player movement without coupling to ECS tilemap components
- Prevents position/collision drift by sharing coordinate system

### Map Format
JSON format (assets/maps/untitled_map.json):
- `'1'` = walkable floor
- `'2'` = solid wall (blocked)
- `'3'` = life support unit (spawns on floor tile)

## Physical Constants

All constants in `src/atmosphere/constants.rs`:

**Gas properties:**
- Universal gas constant: R = 8.314 J/(mol·K)
- Molar masses: M_O₂ = 0.032 kg/mol, M_N₂ = 0.028 kg/mol, M_CO₂ = 0.044 kg/mol
- Earth atmosphere: 101,325 Pa at 21% O₂, 78% N₂, 0.04% CO₂

**Life support system:**
- Room dimensions: 1m × 1m × 2.5m per tile
- Target: 60-second pressurization from vacuum
- Generation rates tuned for 29×14 tile room (1,015 m³)
- Fan radius: 3 tiles with 1.5 m/s tangential swirl speed

**Human physiology:**
- O₂ consumption: 5.54×10⁻⁶ kg/s (~250 mL/min at rest)
- CO₂ production: 6.09×10⁻⁶ kg/s (respiratory quotient 0.8)
- Breathable conditions: 16-50 kPa O₂ partial pressure, <5 kPa CO₂

## Visualization Modes

Press `V` to cycle atmosphere debug overlay:
- **Off**: No visualization
- **Breathability** (default): Neon green = safe air, neon red = dangerous (low O₂ or high CO₂)

Other modes available in code but not exposed to UI:
- Pressure (blue = vacuum, red = 2× Earth pressure)
- Temperature (blue = cold, red = hot)
- Oxygen density (green intensity)
- Nitrogen density (blue intensity)
- CO₂ density (red intensity = danger)
- Velocity field with arrow overlays

## Design Philosophy

This is a **physics simulation testbed** on an accelerated timescale. The goal is realistic atmospheric modeling with:
- Real physical constants (gas laws, diffusion coefficients, human respiration rates)
- Time compression where needed (gas mixing would take days with pure diffusion)
- Physically plausible shortcuts (forced convection fans simulate ventilation systems)
- No arbitrary "gamification" of life support mechanics

## Development Workflow
- **Formatting**: `cargo fmt` before commits
- **Linting**: `cargo clippy --all-targets --all-features`
- **Profiling**: `cargo run --release --features bevy/dynamic_linking`
- **Build checks**: `cargo check` for fast compilation verification

## Known Limitations
- Walls have zero thickness in collision (pure point/AABB detection)
- Player uses point collision instead of bounding box (can clip into walls visually)
- No multi-room pressure differential or breach mechanics yet
- Temperature affects pressure but no heat transfer to/from walls
- Single atmospheric grid (no z-layers or 3D flow)

## Future Considerations
- Doors with pressure seals and airlock cycling
- Hull breach mechanics with rapid decompression
- Multi-room pressure zones with flow between chambers
- Fire simulation (combustion consuming O₂, producing CO₂ and heat)
- HVAC duct networks for explicit air routing
- Hypoxia/hypercapnia player status effects

## Contributing

See `CONTRIBUTING.md` for development workflow, commit conventions, and code standards. When touching atmosphere constants or solver logic, cross-reference the relevant doc in `docs/` and update it so the README stays aligned with the actual implementation.

Key guidelines:
- Keep physics constants in `atmosphere/constants.rs`
- Document units in comments (kg, m, s, Pa, K)
- Test with visualization overlay (V key) to verify gas behavior
- Maintain separation between rendering (tilemap) and collision (collision map)

## License

This project is released under the MIT License. See `LICENSE` (to be added) for full terms.
