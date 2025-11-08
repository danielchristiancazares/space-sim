use super::grid::AtmosphereGrid;
use super::steps::{
    advection::advection_step, diffusion::diffusion_step, pressure_flux::pressure_flux_step,
    pressure_projection::pressure_correction_step,
};
use crate::atmosphere::constants;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::{Res, ResMut, Time};

pub fn simulate_atmosphere(
    mut atmosphere: ResMut<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    time: Res<Time>,
) {
    let mut remaining = time.delta_secs();

    if remaining <= 0.0 {
        return;
    }

    while remaining > 0.0 {
        // For stability, limit timestep based on CFL condition
        let max_velocity = compute_max_velocity(&atmosphere);
        let cfl_dt = if max_velocity > 0.0 {
            constants::CFL_NUMBER * atmosphere.tile_size_physical / max_velocity
        } else {
            remaining
        };

        // Avoid infinitesimal steps that would stall the loop.
        let dt_sim = remaining.min(cfl_dt.max(f32::EPSILON));

        // Operator splitting approach per sub-step:
        advection_step(&mut atmosphere, &collision_map, dt_sim);
        diffusion_step(&mut atmosphere, &collision_map, dt_sim);

        // Update pressure so the subsequent steps see the latest state
        update_grid_pressures(&mut atmosphere);

        pressure_flux_step(&mut atmosphere, &collision_map, dt_sim);

        // Update pressure again after mass transfer
        update_grid_pressures(&mut atmosphere);

        pressure_correction_step(&mut atmosphere, &collision_map, dt_sim);

        remaining -= dt_sim;
    }
}

pub fn update_grid_pressures(atmosphere: &mut AtmosphereGrid) {
    for cell in atmosphere.cells.iter_mut() {
        cell.update_pressure();
    }
}

fn compute_max_velocity(atmosphere: &AtmosphereGrid) -> f32 {
    atmosphere
        .cells
        .iter()
        .map(|cell| (cell.u * cell.u + cell.v * cell.v).sqrt())
        .fold(0.0, f32::max)
}
