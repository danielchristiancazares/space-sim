use crate::atmosphere::constants;
use crate::atmosphere::grid::AtmosphereGrid;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::{debug, warn, Res, ResMut, Resource, Time};

#[derive(Resource, Default)]
pub struct MassTracker {
    pub last_total_mass: f32,
    pub time_since_check: f32,
    pub expected_delta: f32,
}

#[derive(Resource, Default)]
pub struct DivergenceTracker {
    pub time_since_check: f32,
}

#[derive(Resource)]
pub struct MomentumTracker {
    pub last_momentum_x: f32,
    pub last_momentum_y: f32,
    pub time_since_check: f32,
}

impl Default for MomentumTracker {
    fn default() -> Self {
        Self {
            last_momentum_x: 0.0,
            last_momentum_y: 0.0,
            time_since_check: 0.0,
        }
    }
}

pub fn check_mass_conservation(
    atmosphere: Res<AtmosphereGrid>,
    mut tracker: ResMut<MassTracker>,
    time: Res<Time>,
) {
    tracker.time_since_check += time.delta_secs();

    if tracker.time_since_check < constants::MASS_CHECK_INTERVAL {
        return;
    }

    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;
    let total_mass: f32 = atmosphere
        .cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum();

    if tracker.last_total_mass > f32::EPSILON {
        let actual_delta = total_mass - tracker.last_total_mass;
        let error = (actual_delta - tracker.expected_delta).abs();
        let relative_error = if tracker.last_total_mass > f32::EPSILON {
            error / tracker.last_total_mass
        } else {
            0.0
        };

        if relative_error > constants::MASS_TOLERANCE {
            warn!(
                "Mass conservation warning: Total mass changed by {:.6} kg (expected {:.6} kg from sources/sinks). Relative error: {:.2}%",
                actual_delta,
                tracker.expected_delta,
                relative_error * 100.0
            );
        } else {
            debug!(
                "Mass conservation check: Total mass = {:.3} kg (delta: {:.6} kg, expected: {:.6} kg)",
                total_mass, actual_delta, tracker.expected_delta
            );
        }
    }

    tracker.last_total_mass = total_mass;
    tracker.time_since_check = 0.0;
    tracker.expected_delta = 0.0;
}

pub fn monitor_divergence(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut tracker: ResMut<DivergenceTracker>,
    time: Res<Time>,
) {
    tracker.time_since_check += time.delta_secs();

    if tracker.time_since_check < constants::DIVERGENCE_CHECK_INTERVAL {
        return;
    }

    let dx = atmosphere.tile_size_physical;
    let mut max_divergence: f32 = 0.0;
    let mut total_divergence = 0.0;
    let mut cell_count = 0;

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);

            let u_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                atmosphere.cells[atmosphere.index(x + 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                atmosphere.cells[atmosphere.index(x - 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let v_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                atmosphere.cells[atmosphere.index(x, y + 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                atmosphere.cells[atmosphere.index(x, y - 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let du_dx = (u_right - u_left) / (2.0 * dx);
            let dv_dy = (v_up - v_down) / (2.0 * dx);
            let div = du_dx + dv_dy;

            max_divergence = max_divergence.max(div.abs());
            total_divergence += div.abs();
            cell_count += 1;
        }
    }

    let avg_divergence = if cell_count > 0 {
        total_divergence / cell_count as f32
    } else {
        0.0
    };

    debug!(
        "Divergence check: max = {:.6} s⁻¹, avg = {:.6} s⁻¹",
        max_divergence, avg_divergence
    );

    if max_divergence > 1.0 {
        warn!(
            "High velocity divergence: {:.3} s⁻¹ (compressible flow transient; \
             may indicate shock, expansion wave, or CFL violation if sustained)",
            max_divergence
        );
    }

    tracker.time_since_check = 0.0;
}

pub fn check_momentum_conservation(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut tracker: ResMut<MomentumTracker>,
    time: Res<Time>,
) {
    tracker.time_since_check += time.delta_secs();

    if tracker.time_since_check < constants::MASS_CHECK_INTERVAL {
        return;
    }

    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

    let mut total_momentum_x = 0.0;
    let mut total_momentum_y = 0.0;

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells[idx];
            let mass = cell.total_density() * tile_volume;

            total_momentum_x += mass * cell.u;
            total_momentum_y += mass * cell.v;
        }
    }

    if tracker.last_momentum_x.abs() > f32::EPSILON || tracker.last_momentum_y.abs() > f32::EPSILON
    {
        let delta_x = total_momentum_x - tracker.last_momentum_x;
        let delta_y = total_momentum_y - tracker.last_momentum_y;
        let delta_magnitude = (delta_x * delta_x + delta_y * delta_y).sqrt();

        // Momentum can change due to pressure forces and boundary interactions,
        // but large systematic drift may indicate numerical issues
        if delta_magnitude > 10.0 {
            debug!(
                "Momentum change: Δp = ({:.3}, {:.3}) kg·m/s (magnitude: {:.3})",
                delta_x, delta_y, delta_magnitude
            );
        }
    }

    tracker.last_momentum_x = total_momentum_x;
    tracker.last_momentum_y = total_momentum_y;
    tracker.time_since_check = 0.0;
}
