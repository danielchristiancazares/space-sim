use bevy::prelude::*;

use crate::animation::WalkCycle;
use crate::tilemap::TileCollisionMap;

#[derive(Component)]
pub struct Player {
    pub speed: f32,
}

#[derive(Component)]
pub struct Direction {
    pub facing: Vec2,
}

#[derive(Component, Default, PartialEq, Clone, Copy, Debug)]
pub enum MovementState {
    #[default]
    Idle,
    Walking,
    Running,
}

pub struct PlayerPlugin;

impl Plugin for PlayerPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, spawn_player)
            .add_systems(Update, player_movement);
    }
}

fn spawn_player(mut commands: Commands, asset_server: Res<AssetServer>) {
    // Load the south-facing rotation as the default sprite
    let texture = asset_server.load("sprites/characters/rotations/south.png");

    commands.spawn((
        Sprite::from_image(texture),
        Transform {
            translation: Vec3::new(0.0, 0.0, 10.0),
            scale: Vec3::splat(2.0), // Scale up for better visibility
            ..default()
        },
        Player { speed: 150.0 },
        Direction {
            facing: Vec2::new(0.0, -1.0), // Start facing south
        },
        MovementState::Idle,
        WalkCycle::default(),
    ));
}

fn player_movement(
    keyboard_input: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    collision_map: Res<TileCollisionMap>,
    mut query: Query<(&Player, &mut Transform, &mut Direction, &mut MovementState)>,
) {
    for (player, mut transform, mut direction, mut state) in &mut query {
        let mut velocity = Vec2::ZERO;

        // 8-directional WASD input
        if keyboard_input.pressed(KeyCode::KeyW) {
            velocity.y += 1.0;
        }
        if keyboard_input.pressed(KeyCode::KeyS) {
            velocity.y -= 1.0;
        }
        if keyboard_input.pressed(KeyCode::KeyA) {
            velocity.x -= 1.0;
        }
        if keyboard_input.pressed(KeyCode::KeyD) {
            velocity.x += 1.0;
        }

        // Update state and direction
        if velocity.length() > 0.0 {
            // Normalize to prevent faster diagonal movement
            velocity = velocity.normalize();
            direction.facing = velocity;

            // Apply movement
            let delta = velocity * player.speed * time.delta_secs();
            // We treat the player as a point collider for now; expand to AABB once sprites need thickness.
            let mut new_pos = transform.translation.truncate();
            let mut moved = false;

            if delta.x != 0.0 {
                let attempted = Vec2::new(new_pos.x + delta.x, new_pos.y);
                if !collision_map.is_world_blocked(attempted) {
                    new_pos.x = attempted.x;
                    moved = true;
                }
            }

            if delta.y != 0.0 {
                let attempted = Vec2::new(new_pos.x, new_pos.y + delta.y);
                if !collision_map.is_world_blocked(attempted) {
                    new_pos.y = attempted.y;
                    moved = true;
                }
            }

            transform.translation.x = new_pos.x;
            transform.translation.y = new_pos.y;

            *state = if moved {
                MovementState::Walking
            } else {
                MovementState::Idle
            };
        } else {
            *state = MovementState::Idle;
        }
    }
}
