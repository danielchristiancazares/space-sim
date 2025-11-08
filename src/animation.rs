use std::time::Duration;

use bevy::prelude::*;
use bevy::time::TimerMode;

use crate::player::{Direction, MovementState, Player};

const DIRECTION_KEYS: [&str; 8] = [
    "west",
    "south-west",
    "south",
    "south-east",
    "east",
    "north-east",
    "north",
    "north-west",
];

#[derive(Clone)]
pub struct DirectionalFrames {
    pub frames: Vec<Handle<Image>>,
}

#[derive(Resource)]
pub struct WalkAnimation {
    pub speed_fps: f32,
    pub frames_by_direction: Vec<DirectionalFrames>,
}

#[derive(Resource)]
pub struct CharacterSprites {
    pub south: Handle<Image>,
    pub south_west: Handle<Image>,
    pub west: Handle<Image>,
    pub north_west: Handle<Image>,
    pub north: Handle<Image>,
    pub north_east: Handle<Image>,
    pub east: Handle<Image>,
    pub south_east: Handle<Image>,
}

#[derive(Component)]
pub struct WalkCycle {
    timer: Timer,
    frame_index: usize,
    current_direction: usize,
}

impl WalkCycle {
    pub fn new(frame_rate: f32) -> Self {
        Self {
            timer: Timer::from_seconds(1.0 / frame_rate, TimerMode::Repeating),
            frame_index: 0,
            current_direction: 2, // Default facing south.
        }
    }

    pub fn ensure_duration(&mut self, frame_rate: f32) {
        let duration = Duration::from_secs_f32(1.0 / frame_rate);
        if self.timer.duration() != duration {
            self.timer.set_duration(duration);
            self.timer.reset();
        }
    }

    pub fn reset(&mut self) {
        self.frame_index = 0;
        self.timer.reset();
    }
}

impl Default for WalkCycle {
    fn default() -> Self {
        Self::new(10.0)
    }
}

pub struct AnimationPlugin;

impl Plugin for AnimationPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, load_character_sprites)
            .add_systems(Update, animate_character_sprite);
    }
}

fn load_character_sprites(mut commands: Commands, asset_server: Res<AssetServer>) {
    let sprites = CharacterSprites {
        south: asset_server.load("sprites/characters/rotations/south.png"),
        south_west: asset_server.load("sprites/characters/rotations/south-west.png"),
        west: asset_server.load("sprites/characters/rotations/west.png"),
        north_west: asset_server.load("sprites/characters/rotations/north-west.png"),
        north: asset_server.load("sprites/characters/rotations/north.png"),
        north_east: asset_server.load("sprites/characters/rotations/north-east.png"),
        east: asset_server.load("sprites/characters/rotations/east.png"),
        south_east: asset_server.load("sprites/characters/rotations/south-east.png"),
    };

    commands.insert_resource(sprites);

    let base_path = "sprites/characters/animations/walking-8-frames";
    let mut frames_by_direction = Vec::with_capacity(DIRECTION_KEYS.len());

    for dir in DIRECTION_KEYS {
        let mut frames = Vec::new();
        for frame_idx in 0..=5 {
            let path = format!("{base_path}/{dir}/frame_{frame_idx:03}.png");
            frames.push(asset_server.load(path));
        }
        frames_by_direction.push(DirectionalFrames { frames });
    }

    commands.insert_resource(WalkAnimation {
        speed_fps: 10.0,
        frames_by_direction,
    });
}

fn direction_to_index(direction: Vec2) -> usize {
    let angle = direction.y.atan2(direction.x);
    let dir_index =
        ((angle + std::f32::consts::PI) / (std::f32::consts::PI * 2.0) * 8.0).round() as i32;
    (((dir_index % 8) + 8) % 8) as usize
}

fn static_sprite_for_index(sprites: &CharacterSprites, index: usize) -> Handle<Image> {
    match index {
        0 => sprites.west.clone(),
        1 => sprites.south_west.clone(),
        2 => sprites.south.clone(),
        3 => sprites.south_east.clone(),
        4 => sprites.east.clone(),
        5 => sprites.north_east.clone(),
        6 => sprites.north.clone(),
        7 => sprites.north_west.clone(),
        _ => sprites.south.clone(),
    }
}

fn animate_character_sprite(
    time: Res<Time>,
    sprites: Res<CharacterSprites>,
    walk_animation: Res<WalkAnimation>,
    mut query: Query<(&Direction, &MovementState, &mut Sprite, &mut WalkCycle), With<Player>>,
) {
    for (direction, state, mut sprite, mut cycle) in &mut query {
        let dir_index = direction_to_index(direction.facing);
        let Some(frames) = walk_animation.frames_by_direction.get(dir_index) else {
            sprite.image = static_sprite_for_index(&sprites, dir_index);
            continue;
        };

        let is_moving = matches!(*state, MovementState::Walking | MovementState::Running);

        cycle.ensure_duration(walk_animation.speed_fps);

        if cycle.current_direction != dir_index {
            cycle.current_direction = dir_index;
            cycle.reset();
        }

        if is_moving && !frames.frames.is_empty() {
            cycle.timer.tick(time.delta());
            if cycle.timer.finished() {
                cycle.frame_index = (cycle.frame_index + 1) % frames.frames.len();
            }
            let frame = frames.frames[cycle.frame_index].clone();
            sprite.image = frame;
        } else {
            if !frames.frames.is_empty() {
                cycle.frame_index = 0;
            }
            cycle.timer.reset();
            sprite.image = static_sprite_for_index(&sprites, dir_index);
        }
    }
}
