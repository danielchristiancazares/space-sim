use bevy::prelude::*;
use bevy_ecs_tilemap::prelude::*;
use serde::Deserialize;

/// Resource tracking life support tile positions
#[derive(Resource, Default)]
pub struct LifeSupportTiles {
    pub positions: Vec<UVec2>,
}

/// Resource holding wall flags so gameplay systems can run collision queries without
/// depending on rendering components; it mirrors the Tilemap transform to avoid drift.
#[derive(Resource, Clone, Debug)]
pub struct TileCollisionMap {
    pub width: u32,
    pub height: u32,
    pub origin: Vec2,
    pub tile_size: Vec2,
    blocked: Vec<bool>,
}

impl TileCollisionMap {
    fn index(&self, x: u32, y: u32) -> Option<usize> {
        if x < self.width && y < self.height {
            Some((y * self.width + x) as usize)
        } else {
            None
        }
    }

    pub fn is_blocked(&self, x: u32, y: u32) -> bool {
        self.index(x, y)
            .map(|idx| *self.blocked.get(idx).unwrap_or(&false))
            .unwrap_or(true)
    }

    /// Converts a world-space point into the tile coordinates we used to spawn the map.
    pub fn world_to_tile(&self, world_pos: Vec2) -> Option<UVec2> {
        let grid_x = ((world_pos.x - self.origin.x) / self.tile_size.x).floor();
        let grid_y = ((world_pos.y - self.origin.y) / self.tile_size.y).floor();

        if !grid_x.is_finite() || !grid_y.is_finite() {
            return None;
        }

        let x = grid_x as i32;
        let y = grid_y as i32;

        if x < 0 || y < 0 {
            None
        } else {
            let ux = x as u32;
            let uy = y as u32;
            if ux < self.width && uy < self.height {
                Some(UVec2::new(ux, uy))
            } else {
                None
            }
        }
    }

    /// Returns true when the point is outside the map or lands on a wall tile.
    pub fn is_world_blocked(&self, world_pos: Vec2) -> bool {
        self.world_to_tile(world_pos)
            .map(|tile| self.is_blocked(tile.x, tile.y))
            .unwrap_or(true)
    }

    /// Create a test collision map with all tiles unblocked (for testing).
    #[cfg(test)]
    pub fn test_map(width: u32, height: u32) -> Self {
        Self {
            width,
            height,
            origin: Vec2::ZERO,
            tile_size: Vec2::new(1.0, 1.0),
            blocked: vec![false; (width * height) as usize],
        }
    }

    /// Create an empty collision map with specified dimensions and transform.
    #[allow(dead_code)]
    pub fn empty(width: u32, height: u32, tile_size: Vec2, origin: Vec2) -> Self {
        Self {
            width,
            height,
            origin,
            tile_size,
            blocked: vec![false; (width * height) as usize],
        }
    }
}

const MAP_JSON: &str = include_str!("../assets/maps/untitled_map.json");

#[derive(Deserialize)]
struct SimpleMap {
    width: u32,
    height: u32,
    rows: Vec<String>,
}

pub struct TilemapPlugin;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub struct TilemapInitSet;

impl Plugin for TilemapPlugin {
    fn build(&self, app: &mut App) {
        app.add_plugins(bevy_ecs_tilemap::TilemapPlugin)
            .add_systems(Startup, create_tilemap.in_set(TilemapInitSet));
    }
}

fn create_tilemap(mut commands: Commands, asset_server: Res<AssetServer>) {
    let mut life_support_positions = Vec::new();
    let map_data: SimpleMap =
        serde_json::from_str(MAP_JSON).expect("untitled map JSON should deserialize");

    // Tile indices for different tile types:
    const FLOOR_TILE_INDEX: u32 = 0;
    const WALL_TILE_INDEX: u32 = 15;
    const LIFE_SUPPORT_TILE_INDEX: u32 = 0; // Will use separate sprite

    let map_size = TilemapSize {
        x: map_data.width,
        y: map_data.height,
    };

    assert_eq!(
        map_data.rows.len() as u32,
        map_size.y,
        "map height should match row count"
    );

    // Using compile-time inclusion keeps startup deterministic without adding a runtime asset loader.
    // Load the 32x32 Wang tileset (4x4 grid = 16 tiles)
    let tileset_texture: Handle<Image> = asset_server.load("sprites/tiles/metal_floor_tileset.png");

    let tile_size = TilemapTileSize { x: 32.0, y: 32.0 };
    let grid_size = TilemapGridSize { x: 32.0, y: 32.0 };
    let map_type = TilemapType::Square;

    let tilemap_entity = commands.spawn_empty().id();

    let mut tile_storage = TileStorage::empty(map_size);
    let mut blocked_tiles = Vec::with_capacity((map_size.x * map_size.y) as usize); // Keep collision data alongside tiles.

    for (y, row) in map_data.rows.iter().enumerate() {
        assert_eq!(
            row.len() as u32,
            map_size.x,
            "row {y} width should match map width"
        );

        let row_bytes = row.as_bytes();

        for (x, cell) in row_bytes.iter().enumerate() {
            let (tile_index, is_wall, is_life_support) = match cell {
                b'1' => (FLOOR_TILE_INDEX, false, false),
                b'2' => (WALL_TILE_INDEX, true, false),
                b'3' => (LIFE_SUPPORT_TILE_INDEX, false, true), // Life support on floor
                _ => (FLOOR_TILE_INDEX, false, false), // Unrecognized cells default to floor
            };

            blocked_tiles.push(is_wall);

            if is_life_support {
                life_support_positions.push(UVec2::new(x as u32, y as u32));
            }

            // Walls all use the solid tile for now; if we add corner sprites later, this match can branch on adjacency.
            let tile_pos = TilePos {
                x: x as u32,
                y: y as u32,
            };
            let tile_entity = commands
                .spawn(TileBundle {
                    position: tile_pos,
                    tilemap_id: TilemapId(tilemap_entity),
                    texture_index: TileTextureIndex(tile_index),
                    ..Default::default()
                })
                .id();

            tile_storage.set(&tile_pos, tile_entity);
        }
    }

    // Center the tilemap at the origin
    // Calculate offset: -(map_size * tile_size) / 2
    let offset_x = -(map_size.x as f32 * tile_size.x) / 2.0;
    let offset_y = -(map_size.y as f32 * tile_size.y) / 2.0;

    commands.entity(tilemap_entity).insert((TilemapBundle {
        grid_size,
        map_type,
        size: map_size,
        storage: tile_storage,
        texture: TilemapTexture::Vector(vec![tileset_texture]),
        tile_size,
        transform: Transform::from_xyz(offset_x, offset_y, 0.0),
        ..Default::default()
    },));

    commands.insert_resource(TileCollisionMap {
        width: map_size.x,
        height: map_size.y,
        origin: Vec2::new(offset_x, offset_y),
        tile_size: Vec2::new(tile_size.x, tile_size.y),
        blocked: blocked_tiles,
    });

    commands.insert_resource(LifeSupportTiles {
        positions: life_support_positions.clone(),
    });

    // Spawn life support sprites
    let life_support_texture: Handle<Image> = asset_server.load("sprites/life_support.png");

    for pos in life_support_positions {
        let world_x = offset_x + (pos.x as f32 + 0.5) * tile_size.x;
        let world_y = offset_y + (pos.y as f32 + 0.5) * tile_size.y;

        commands.spawn((
            Sprite::from_image(life_support_texture.clone()),
            Transform::from_xyz(world_x, world_y, 1.0), // Z=1 to be above floor
        ));
    }
}
