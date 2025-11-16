use crate::player::Player;
use bevy::prelude::*;

#[derive(Component)]
pub struct MainCamera {
    pub smoothness: f32,
}

pub struct CameraPlugin;

impl Plugin for CameraPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, spawn_camera)
            .add_systems(Update, camera_follow_player);
    }
}

fn spawn_camera(mut commands: Commands) {
    commands.spawn((Camera2d, MainCamera { smoothness: 5.0 }));
}

fn camera_follow_player(
    time: Res<Time>,
    player_query: Query<&Transform, With<Player>>,
    mut camera_query: Query<(&mut Transform, &MainCamera), (Without<Player>, With<Camera2d>)>,
) {
    if let Ok(player_transform) = player_query.single() {
        for (mut camera_transform, camera) in &mut camera_query {
            // Smooth interpolation towards player position
            let target = player_transform.translation;
            let current = camera_transform.translation;

            let lerp_factor = camera.smoothness * time.delta_secs();
            let new_pos = current.lerp(target, lerp_factor.min(1.0));

            camera_transform.translation =
                Vec3::new(new_pos.x, new_pos.y, camera_transform.translation.z);
        }
    }
}
