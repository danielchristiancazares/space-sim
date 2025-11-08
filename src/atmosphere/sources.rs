use super::grid::AtmosphereGrid;
use super::monitoring::MassTracker;
use crate::atmosphere::constants;
use crate::player::Player;
use crate::tilemap::{LifeSupportTiles, TileCollisionMap};
use bevy::prelude::*;

pub fn life_support_generation(
    mut atmosphere: ResMut<AtmosphereGrid>,
    life_support: Res<LifeSupportTiles>,
    mut tracker: ResMut<MassTracker>,
    time: Res<Time>,
) {
    let dt = time.delta_secs();

    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

    let mut total_generated = 0.0;

    for pos in &life_support.positions {
        let idx = atmosphere.index(pos.x, pos.y);
        let cell = &mut atmosphere.cells[idx];

        let o2_generated = constants::LIFE_SUPPORT_O2_RATE * dt;
        let n2_generated = constants::LIFE_SUPPORT_N2_RATE * dt;

        let delta_rho_o2 = o2_generated / tile_volume;
        let delta_rho_n2 = n2_generated / tile_volume;

        cell.rho_o2 += delta_rho_o2;
        cell.rho_n2 += delta_rho_n2;

        total_generated += o2_generated + n2_generated;

        if cell.temperature < constants::ROOM_TEMP {
            cell.temperature = constants::ROOM_TEMP;
        }
    }

    tracker.expected_delta += total_generated;
}

pub fn life_support_mixing(
    mut atmosphere: ResMut<AtmosphereGrid>,
    life_support: Res<LifeSupportTiles>,
    collision_map: Res<TileCollisionMap>,
) {
    if life_support.positions.is_empty() {
        return;
    }

    let radius = constants::FAN_RADIUS_TILES;
    let radius_f32 = radius as f32;

    for pos in &life_support.positions {
        let center_x = pos.x as i32;
        let center_y = pos.y as i32;

        for dy in -radius..=radius {
            for dx in -radius..=radius {
                if dx == 0 && dy == 0 {
                    continue;
                }

                let x = center_x + dx;
                let y = center_y + dy;

                if !atmosphere.in_bounds(x, y) || collision_map.is_blocked(x as u32, y as u32) {
                    continue;
                }

                let Some(cell) = atmosphere.get_mut(x, y) else {
                    continue;
                };

                let offset = Vec2::new(dx as f32, dy as f32);
                let distance = offset.length();
                if distance > radius_f32 {
                    continue;
                }

                let tangent = Vec2::new(-offset.y, offset.x);
                let tangent_len = tangent.length();
                if tangent_len <= f32::EPSILON {
                    continue;
                }
                let swirl_dir = tangent / tangent_len;

                let falloff = (1.0 - distance / (radius_f32 + 1.0)).clamp(0.0, 1.0);
                if falloff <= 0.0 {
                    continue;
                }

                let target_velocity = swirl_dir * constants::FAN_SWIRL_SPEED;
                let blend = constants::FAN_BLEND * falloff;

                cell.u += (target_velocity.x - cell.u) * blend;
                cell.v += (target_velocity.y - cell.v) * blend;
            }
        }
    }
}

pub fn player_respiration(
    mut atmosphere: ResMut<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut tracker: ResMut<MassTracker>,
    player_query: Query<&Transform, With<Player>>,
    time: Res<Time>,
) {
    let Ok(player_transform) = player_query.get_single() else {
        return;
    };

    let player_pos = player_transform.translation.truncate();

    let Some(tile_pos) = collision_map.world_to_tile(player_pos) else {
        return;
    };

    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

    let idx = atmosphere.index(tile_pos.x, tile_pos.y);
    let cell = &mut atmosphere.cells[idx];

    let dt = time.delta_secs();
    let o2_consumed = constants::O2_CONSUMPTION_RATE * dt;
    let co2_produced = constants::CO2_PRODUCTION_RATE * dt;

    let delta_rho_o2 = -o2_consumed / tile_volume;
    let delta_rho_co2 = co2_produced / tile_volume;

    cell.rho_o2 = (cell.rho_o2 + delta_rho_o2).max(0.0);
    cell.rho_co2 += delta_rho_co2;
    cell.update_pressure();

    let net_mass_change = co2_produced - o2_consumed;
    tracker.expected_delta += net_mass_change;

    let o2_partial_pressure = (cell.rho_o2 / constants::M_O2) * constants::R * cell.temperature;
    if cell.pressure > 0.0 && o2_partial_pressure < 16000.0 {
        warn!(
            "Low oxygen at player location! O2 partial pressure: {:.0} Pa ({:.1}%)",
            o2_partial_pressure,
            100.0 * o2_partial_pressure / cell.pressure
        );
    }

    let co2_partial_pressure = (cell.rho_co2 / constants::M_CO2) * constants::R * cell.temperature;
    if cell.pressure > 0.0 && co2_partial_pressure > 5000.0 {
        warn!(
            "High CO2 at player location! CO2 partial pressure: {:.0} Pa ({:.1}%)",
            co2_partial_pressure,
            100.0 * co2_partial_pressure / cell.pressure
        );
    }
}
