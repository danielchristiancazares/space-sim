use bevy::log::LogPlugin;
use bevy::prelude::*;

mod animation;
mod atmosphere;
mod camera;
mod player;
mod tilemap;

use animation::AnimationPlugin;
use atmosphere::AtmospherePlugin;
use camera::CameraPlugin;
use player::PlayerPlugin;
use tilemap::TilemapPlugin;

fn main() {
    App::new()
        .add_plugins(
            DefaultPlugins
                .set(ImagePlugin::default_nearest())
                .set(WindowPlugin {
                    primary_window: Some(Window {
                        title: "Space Station".to_string(),
                        resolution: (1280, 720).into(),
                        ..default()
                    }),
                    ..default()
                })
                .set(LogPlugin {
                    // Research mode: Clean output without timestamps
                    filter: "info,wgpu_core=warn,wgpu_hal=warn,bevy_render::renderer=warn,bevy_winit=warn".into(),
                    level: bevy::log::Level::INFO,
                    ..default()
                }),
        )
        .add_plugins((
            TilemapPlugin,
            PlayerPlugin,
            CameraPlugin,
            AnimationPlugin,
            AtmospherePlugin,
        ))
        .run();
}
