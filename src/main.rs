use bevy::prelude::*;

mod animation;
mod atmosphere;
mod camera;
mod debug;
mod player;
mod tilemap;

use animation::AnimationPlugin;
use atmosphere::AtmospherePlugin;
use camera::CameraPlugin;
use debug::PressureLoggerPlugin;
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
                        resolution: (1280.0, 720.0).into(),
                        ..default()
                    }),
                    ..default()
                }),
        )
        .add_plugins((
            TilemapPlugin,
            PlayerPlugin,
            CameraPlugin,
            AnimationPlugin,
            AtmospherePlugin,
            PressureLoggerPlugin,
        ))
        .run();
}
