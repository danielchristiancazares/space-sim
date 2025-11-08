# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

A 2D space station game built with Bevy 0.15 game engine featuring 8-directional player movement, smooth camera following, and tile-based collision detection with `bevy_ecs_tilemap`.

## Git Workflow

### Branch Naming Convention

**IMPORTANT**: Use conventional topic branch names, NOT `claude/*` prefixes.

When creating branches, follow these conventions:
- `fix/description` - For bug fixes (e.g., `fix/momentum-physics`)
- `feature/description` - For new features (e.g., `feature/life-support`)
- `refactor/description` - For code refactoring (e.g., `refactor/atmosphere-grid`)
- `docs/description` - For documentation updates (e.g., `docs/physics-constants`)
- `test/description` - For adding tests (e.g., `test/conservation-checks`)

### Examples

```bash
# Good ✓
git checkout -b fix/pressure-correction
git checkout -b feature/airlock-system
git checkout -b refactor/collision-detection

# Bad ✗
git checkout -b claude/fix-physics
git checkout -b my-changes
git checkout -b update
```

## Build Commands

```bash
# Development build (fast iteration with dependency optimizations)
cargo run

# Release build (optimized with thin LTO)
cargo run --release

# Check compilation without running
cargo check

# Build without running
cargo build
cargo build --release

# Clean build artifacts
cargo clean
```

## Architecture Overview

### Plugin-Based System

The game uses Bevy's plugin architecture with four main plugins registered in `src/main.rs:27-32`:

1. **TilemapPlugin** - Handles tile rendering and collision map initialization
2. **PlayerPlugin** - Player spawning and movement logic
3. **CameraPlugin** - Camera spawning and smooth following
4. **AnimationPlugin** - Character sprite rotation based on direction

### Key Architectural Patterns

**Collision System**: The game separates rendering from gameplay logic through `TileCollisionMap` (src/tilemap.rs:7-14):
- A resource containing a grid of blocked/walkable flags
- Mirrors the visual tilemap's transform (origin, tile_size) to convert world positions to tile coordinates
- Queried by player movement system (src/player.rs:91-99) without depending on tilemap ECS components
- Prevents drift between visual and collision representations

**Map Loading**: Maps are loaded from JSON at compile-time using `include_str!` (src/tilemap.rs:64):
- Format: `{"width": 31, "height": 16, "rows": ["222...1112", ...]}`
- Character codes: `'1'` = floor (walkable), `'2'` = wall (blocked)
- Located in `assets/maps/untitled_map.json`
- Parsed during `create_tilemap` system and used to initialize both visual tiles and collision map

**Movement System**: The player movement (src/player.rs:54-117) uses per-axis collision checking:
- Attempts X movement first, only applying if not blocked
- Then attempts Y movement independently
- Allows sliding along walls when moving diagonally into obstacles
- Uses normalized velocity to prevent diagonal speed advantage

**Direction and Animation**: Character sprite rotation (src/animation.rs:43-69):
- Direction component stores facing vector on player entity
- Animation system calculates angle from facing vector using `atan2`
- Converts angle to 8-directional index (0-7) representing compass directions
- Sprite changes are tracked via `Changed<Direction>` query filter for efficiency
- Sprite handles loaded once at startup into `CharacterSprites` resource

**Camera System**: Smooth interpolation (src/camera.rs:25-42):
- Queries player Transform and camera Transform separately (Without<Player> filter prevents overlap)
- Uses `lerp` with time-scaled factor for framerate-independent smoothing
- Smoothness parameter (default 5.0) controls follow speed

## Asset Structure

```
assets/
├── maps/
│   └── untitled_map.json          # Tile layout with collision data
└── sprites/
    ├── characters/
    │   └── rotations/              # 8 directional sprites (north.png, south-east.png, etc.)
    └── tiles/
        └── metal_floor_tileset.png # 32x32 Wang tileset (4x4 grid = 16 tiles)
```

### Adding New Maps

1. Create JSON file in `assets/maps/` with format:
   ```json
   {"width": 31, "height": 16, "rows": ["2111...112", ...]}
   ```
2. Update `MAP_JSON` constant in src/tilemap.rs:64 to point to new file
3. Map must be rectangular (all rows same length)
4. Use `'1'` for walkable floor, `'2'` for solid walls

### Adding New Tile Types

Currently the tilemap uses hardcoded indices (src/tilemap.rs:90-91):
- `FLOOR_TILE_INDEX = 0` for walkable tiles
- `WALL_TILE_INDEX = 15` for solid walls

To add more tile types:
1. Add new constants for tile indices
2. Update the cell matching logic in src/tilemap.rs:128-132
3. Add new character codes to the map JSON format
4. Update collision logic if new tiles have different blocking behavior

## Important Configuration

**Window Settings**: Configured in src/main.rs:18-24 with nearest-neighbor filtering for pixel-perfect rendering

**Bevy Version**: Using Bevy 0.15 and bevy_ecs_tilemap 0.15 - both must stay in sync

**Build Optimizations**: Cargo.toml:13-23 sets aggressive dependency optimization in dev builds for faster iteration while keeping main code compilation fast
- The "game" is less a game and more of a physics simulation test on an accelerated timescale.