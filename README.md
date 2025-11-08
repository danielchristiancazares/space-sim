# Space Station Atmospheric Simulation

A 2D fluid dynamics physics simulation built with Bevy 0.15, featuring realistic atmospheric modeling with compressible Navier-Stokes fluid dynamics. The project simulates O₂/N₂/CO₂ gas diffusion, pressure dynamics, player respiration, and life support systems with forced convection on an accelerated timescale. This repository is a simulation testbed rather than a traditional game.

## Feature Overview

### Atmospheric Physics
- **Compressible Navier-Stokes simulation** with operator splitting (advection, diffusion, pressure correction)
- **Multi-gas tracking**: Oxygen (O₂), Nitrogen (N₂), and Carbon Dioxide (CO₂) with separate density fields
- **Gas diffusion** with Fick's law and realistic diffusion coefficients
- **Pressure dynamics** computed from ideal gas law: p = (ρ_O₂/M_O₂ + ρ_N₂/M_N₂ + ρ_CO₂/M_CO₂) × R × T
- **Temperature simulation** with thermal diffusion
- **CFL-stable timestep** sub-stepping for numerical stability
- **Life support systems** with forced convection fans creating swirl patterns for gas mixing
- **Player respiration** consuming O₂ and producing CO₂ at realistic rates

### Rendering Engine Systems
- **8-directional character movement** with normalized diagonal speed and per-axis collision
- **Smooth camera following** with configurable interpolation
- **Tile-based collision detection** using decoupled collision map
- **Character sprite rotation** for 8 compass directions
- **Real-time visualization** (press V): Pressure, temperature, O₂, N₂, CO₂, breathability, and velocity fields

## Technology Stack
- [Bevy v0.15](https://docs.rs/bevy/0.15.0/bevy/) – ECS game engine providing rendering, input, and scheduling
- [bevy_ecs_tilemap v0.15](https://docs.rs/bevy_ecs_tilemap/0.15.0/bevy_ecs_tilemap/) – GPU-accelerated tilemap renderer
- Rust – Safe systems language for deterministic physics simulation

## Getting Started

### Prerequisites
1. Install Rust stable toolchain via `rustup`
2. Ensure `cargo` is on your PATH

### Installation
Clone the repository - no additional dependencies required (Bevy bundles rendering backends).

## Running the Simulation

Development build (faster compile, optimized dependencies):
```bash
cargo run
```

Release build (optimized simulation performance):
```bash
cargo run --release
```

## Controls
- `W` / `A` / `S` / `D` – Move (8-directional with diagonals)
- `V` – Toggle atmosphere visualization (cycles breathability overlay on/off)

## Project Structure
```
space-sim/
├── assets/
│   ├── maps/
│   │   └── untitled_map.json       # Tile layout with '3' marking life support units
│   └── sprites/
│       ├── characters/
│       │   └── rotations/          # 8-directional character sprites
│       ├── tiles/
│       │   └── metal_floor_tileset.png  # 32×32 Wang tileset
│       └── life_support.png
├── src/
│   ├── animation.rs                # Character sprite rotation system
│   ├── atmosphere/                 # Atmospheric simulation modules
│   │   ├── constants.rs            # Physical/tuning constants
│   │   ├── grid.rs                 # Atmosphere grid resource & cell definitions
│   │   ├── monitoring.rs           # Mass/divergence trackers
│   │   ├── sources.rs              # Life support & respiration systems
│   │   ├── steps/                  # CFD operator implementations
│   │   │   ├── advection.rs
│   │   │   ├── diffusion.rs
│   │   │   ├── pressure_flux.rs
│   │   │   └── pressure_projection.rs
│   │   ├── debug.rs                # Visualization helpers
│   │   ├── plugin.rs               # Bevy plugin + schedule wiring
│   │   └── simulation.rs           # Frame orchestration
│   ├── camera.rs                   # Smooth camera follow
│   ├── debug.rs                    # Pressure logger
│   ├── main.rs                     # Plugin registration
│   ├── player.rs                   # Movement and collision
│   └── tilemap.rs                  # Map loading and collision grid
├── Cargo.toml
└── README.md
```

## Architecture Overview

### Atmospheric Simulation Pipeline
The atmosphere plugin orchestrates a four-step CFD loop each frame:

1. **Life support generation** (`sources.rs`): Adds O₂/N₂ at life support tiles based on `LIFE_SUPPORT_O2_RATE` (4.7 kg/s) and `LIFE_SUPPORT_N2_RATE` (15.6 kg/s)

2. **Forced convection** (`sources.rs`): Fans create 1.5 m/s tangential swirl within a 3-tile radius to distribute gases faster than diffusion alone

3. **Player respiration** (`sources.rs`): Consumes O₂ at 5.54×10⁻⁶ kg/s and produces CO₂ at 6.09×10⁻⁶ kg/s (human resting rate)

4. **Pressure update** (`simulation.rs`): Recomputes pressure from densities using the ideal gas law

5. **Fluid dynamics simulation** (`steps/`):
   - **Advection** (`steps/advection.rs`): MacCormack transport of density, momentum, and temperature
   - **Diffusion** (`steps/diffusion.rs`): Viscous and scalar diffusion
   - **Pressure-driven flux** (`steps/pressure_flux.rs`): Mass exchange along pressure gradients with momentum transfer
   - **Pressure projection** (`steps/pressure_projection.rs`): Enforces near-zero velocity divergence

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

## Contribution Guidelines
- Keep physics constants in `atmosphere/constants.rs`
- Document units in comments (kg, m, s, Pa, K)
- Test with visualization overlay (V key) to verify gas behavior
- Maintain separation between rendering (tilemap) and collision (collision map)

## License

This project is released under the MIT License. See `LICENSE` (to be added) for full terms.
