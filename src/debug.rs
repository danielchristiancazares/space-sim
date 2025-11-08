use bevy::prelude::*;
use std::fs::{create_dir_all, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use crate::atmosphere::{AtmosphereGrid, AtmosphereSimSet};
use crate::player::Player;
use crate::tilemap::TileCollisionMap;

/// Plugin that periodically samples the atmospheric pressure at the player's tile and appends it
/// to a CSV for offline analysis.
pub struct PressureLoggerPlugin;

#[derive(Resource)]
struct PressureLogConfig {
    path: PathBuf,
}

#[derive(Resource)]
struct PressureLogState {
    timer: Timer,
}

impl Plugin for PressureLoggerPlugin {
    fn build(&self, app: &mut App) {
        app.insert_resource(PressureLogConfig {
            path: PathBuf::from("debug/pressure.csv"),
        })
        .insert_resource(PressureLogState {
            timer: Timer::from_seconds(1.0, TimerMode::Repeating),
        })
        .add_systems(Update, sample_pressure.after(AtmosphereSimSet::Main));
    }
}

fn sample_pressure(
    config: Res<PressureLogConfig>,
    player_query: Query<&Transform, With<Player>>,
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut state: ResMut<PressureLogState>,
    time: Res<Time>,
) {
    if !state.timer.tick(time.delta()).just_finished() {
        return;
    }

    let Ok(player_transform) = player_query.get_single() else {
        return;
    };

    let Some(tile_pos) = collision_map.world_to_tile(player_transform.translation.truncate())
    else {
        return;
    };

    let idx = atmosphere.index(tile_pos.x, tile_pos.y);
    let pressure = atmosphere.cells[idx].pressure;

    if let Some(parent) = config.path.parent() {
        if let Err(err) = create_dir_all(parent) {
            error!(
                "Failed to create pressure log directory {:?}: {err}",
                parent
            );
            return;
        }
    }

    let mut file = match OpenOptions::new()
        .create(true)
        .append(true)
        .open(&config.path)
    {
        Ok(file) => file,
        Err(err) => {
            error!("Failed to open pressure log {:?}: {err}", config.path);
            return;
        }
    };

    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or(Duration::ZERO);
    if let Err(err) = writeln!(file, "{:.3},{}", timestamp.as_secs_f64(), pressure) {
        error!("Failed to write pressure log entry: {err}");
    }
}
