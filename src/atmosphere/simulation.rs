use super::grid::AtmosphereGrid;
use super::observability::{CsvExporter, SimulationDiagnostics};
use super::steps::{
    advection::advection_step, compression::compression_heating_step, diffusion::diffusion_step,
};
use crate::atmosphere::constants;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::{Res, ResMut, Time};

fn compute_total_mass(atmosphere: &AtmosphereGrid) -> f32 {
    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;
    atmosphere
        .cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum()
}

fn ensure_state_finite(stage: &str, atmosphere: &AtmosphereGrid) {
    for (idx, cell) in atmosphere.cells.iter().enumerate() {
        if !cell.u.is_finite()
            || !cell.v.is_finite()
            || !cell.temperature.is_finite()
            || !cell.rho_o2.is_finite()
            || !cell.rho_n2.is_finite()
            || !cell.rho_co2.is_finite()
        {
            let x = idx % atmosphere.width as usize;
            let y = idx / atmosphere.width as usize;
            panic!(
                "{stage}: Non-finite cell at ({}, {}): u={}, v={}, T={}, rho_o2={}, rho_n2={}, rho_co2={}",
                x, y, cell.u, cell.v, cell.temperature, cell.rho_o2, cell.rho_n2, cell.rho_co2
            );
        }
    }
}

pub fn simulate_atmosphere(
    mut atmosphere: ResMut<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut diagnostics: ResMut<SimulationDiagnostics>,
    mut csv_exporter: ResMut<CsvExporter>,
    time: Res<Time>,
) {
    // FixedUpdate provides consistent dt, but we still check CFL for safety
    let dt = time.delta_secs();

    if dt <= 0.0 {
        return;
    }

    // Acoustic CFL condition for compressible flow
    // Wave speed = |u| + c where c = √(γRT) is sound speed
    // CFL constraint: (|u| + c)Δt/Δx ≤ CFL_NUMBER ensures both advective
    // and acoustic waves are properly resolved.
    let max_wave_speed = compute_max_wave_speed(&atmosphere);

    let cfl_dt = if max_wave_speed > 0.0 {
        constants::CFL_NUMBER * atmosphere.tile_size_physical / max_wave_speed
    } else {
        dt
    };

    // Research-grade accuracy: Follow CFL condition precisely (no artificial cap)
    // This may result in many substeps for high wave speeds, but ensures stability
    let num_substeps = ((dt / cfl_dt).ceil() as usize).max(1);
    let dt_sim = dt / num_substeps as f32;

    bevy::log::trace!(
        "CFL: max_wave_speed={:.1} m/s, cfl_dt={:.4} ms, substeps={}, dt_sim={:.4} ms",
        max_wave_speed,
        cfl_dt * 1000.0,
        num_substeps,
        dt_sim * 1000.0
    );

    for substep in 0..num_substeps {
        let mass_before = compute_total_mass(&atmosphere);

        // Operator splitting approach per sub-step:
        advection_step(&mut atmosphere, &collision_map, dt_sim);
        ensure_state_finite("after advection", &atmosphere);

        let mass_after_advection = compute_total_mass(&atmosphere);
        // Mass rescaling removed for research accuracy - conservation must be inherent

        diffusion_step(&mut atmosphere, &collision_map, dt_sim);
        ensure_state_finite("after diffusion", &atmosphere);

        let mass_after_diffusion = compute_total_mass(&atmosphere);
        // Mass rescaling removed for research accuracy - conservation must be inherent

        // Compression/expansion: Applies -(γ-1)T∇·u to temperature based on local divergence.
        // Density compression (∂ρ/∂t = -ρ∇·u) remains embedded in the flux-conservative advection step.
        compression_heating_step(&mut atmosphere, &collision_map, dt_sim);
        ensure_state_finite("after compression", &atmosphere);

        let mass_after_compression = compute_total_mass(&atmosphere);

        // Update pressure from new ρ and T (EOS consistency)
        update_grid_pressures(&mut atmosphere);

        if substep == 0 || substep == num_substeps - 1 {
            bevy::log::debug!(
                "Substep {}/{}: dt={:.4}, mass: before={:.6}, adv={:.6} (Δ={:.6}), diff={:.6} (Δ={:.6}), \
                 comp={:.6} (Δ={:.6})",
                substep + 1,
                num_substeps,
                dt_sim,
                mass_before,
                mass_after_advection,
                mass_after_advection - mass_before,
                mass_after_diffusion,
                mass_after_diffusion - mass_after_advection,
                mass_after_compression,
                mass_after_compression - mass_after_diffusion,
            );
        }
    }

    // Update diagnostics after all substeps
    diagnostics.update(&atmosphere, &collision_map, dt);

    // Export to CSV at configured intervals
    csv_exporter.update(&diagnostics, dt);
}

pub fn update_grid_pressures(atmosphere: &mut AtmosphereGrid) {
    for cell in atmosphere.cells.iter_mut() {
        cell.update_pressure();
    }
}

/// Compute maximum acoustic wave speed: |u| + c
///
/// For compressible flow, waves propagate at the characteristic speed:
/// - Advection: velocity |u|
/// - Acoustic waves: sound speed c = √(γRT/M)
///
/// Total wave speed = |u| + c ensures CFL condition captures both phenomena.
///
/// Returns: Maximum wave speed [m/s]
fn compute_max_wave_speed(atmosphere: &AtmosphereGrid) -> f32 {
    atmosphere
        .cells
        .iter()
        .map(|cell| {
            let velocity = (cell.u * cell.u + cell.v * cell.v).sqrt();

            // Prevent sqrt(0) at vacuum/CMB temperatures
            let temperature = cell.temperature.max(constants::T_CMB);

            // Sound speed: c = √(γ × R_specific × T)
            // For air at 293K: c ≈ √(1.4 × 287 × 293) ≈ 343 m/s
            let sound_speed = (constants::GAMMA * constants::R_SPECIFIC_AIR * temperature).sqrt();

            // Total characteristic wave speed
            velocity + sound_speed
        })
        .fold(0.0, f32::max)
}

// Periodic violation logging and detailed stats are provided by observability::monitors
