//! Compression/expansion heating for temperature.
//!
//! Implements the energy equation compression term:
//! ∂T/∂t = -(γ-1)T∇·u
//!
//! Note: Density compression (∂ρ/∂t = -ρ∇·u) is handled by flux-conservative
//! advection in the advection module, not here.
//!
//! Physical interpretation:
//! - When gas compresses (∇·u < 0): Temperature increases (compression heating)
//! - When gas expands (∇·u > 0): Temperature decreases (expansion cooling)
//!
//! This is the thermodynamic coupling that makes adiabatic processes realistic.

use crate::atmosphere::constants;
use crate::atmosphere::grid::AtmosphereGrid;
use crate::tilemap::TileCollisionMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Apply compression/expansion heating to temperature field.
///
/// Implements:
///   ∂T/∂t = -(γ-1)T∇·u         (energy equation - compression heating)
///
/// Where:
/// - T = temperature [K]
/// - γ = Cp/Cv ≈ 1.4 (ratio of specific heats for air)
/// - ∇·u = velocity divergence [s⁻¹]
///
/// Physical examples:
/// - Gas expanding into vacuum: ∇·u > 0 → T↓ (expansion cooling)
/// - Gas compressed by piston: ∇·u < 0 → T↑ (compression heating)
/// - Spray can: Rapid expansion → feels cold (temperature drops)
pub fn compression_heating_step(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;
    let width = atmosphere.width as usize;

    // Counters for hot-spot brake diagnostics (thread-safe)
    let hot_zone_entries = AtomicUsize::new(0);
    let phys_max_caps = AtomicUsize::new(0);

    // Buffer current state for divergence calculation
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    let cells_buffer = &atmosphere.cells_buffer;
    let cells = &mut atmosphere.cells;

    cells
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(y, row)| {
            let y_u32 = y as u32;

            for x in 0..width {
                let x_u32 = x as u32;

                if collision_map.is_blocked(x_u32, y_u32) {
                    continue;
                }

                let idx = y * width + x;
                let cell = &cells_buffer[idx];

                // Compute velocity divergence: ∇·u = ∂u/∂x + ∂v/∂y
                let div_u_raw = compute_divergence(
                    cells_buffer,
                    atmosphere.width,
                    atmosphere.height,
                    collision_map,
                    x_u32,
                    y_u32,
                    dx,
                );

                // No divergence clamp - CFL enforcement ensures stability
                let div_u = div_u_raw;

                // Guard against extreme divergence that would indicate CFL violation
                debug_assert!(
                    div_u.abs() < 1000.0,
                    "Extreme divergence at ({}, {}): div_u = {:.2} s⁻¹ \
                     (indicates CFL violation or numerical instability)",
                    x_u32,
                    y_u32,
                    div_u
                );

                let t_old = cell.temperature.max(constants::T_CMB);
                let gamma_minus_1 = constants::GAMMA - 1.0;

                // ============================================================
                // Temperature compression: ∂T/∂t = -(γ-1)T∇·u (energy equation)
                // ============================================================
                // EXACT EXPONENTIAL UPDATE with asymmetric clamping and hot-spot brake:
                //   T(t+dt) = T(t) × exp[-(γ-1) × ∇·u × dt]
                //
                // Work in log-space for better numerical control:
                //   dlogT = -(γ-1) × ∇·u × dt

                // Step 1: Compute raw log-change
                let dlog_t_raw = -gamma_minus_1 * div_u * dt;

                // Step 2: Per-second fractional limits, scaled by dt
                // This ensures rate limits are independent of substep size
                let max_log_pos = (1.0 + constants::MAX_FRAC_T_POS_PER_S * dt).ln();
                let max_log_neg = (1.0 - constants::MAX_FRAC_T_NEG_PER_S * dt).ln();

                let mut dlog_t = dlog_t_raw.clamp(max_log_neg, max_log_pos);

                // Step 3: Zone-based hot-spot brake (applied only to heating)
                if dlog_t > 0.0 {
                    if t_old >= constants::T_PHYS_MAX {
                        // Zone 3: No heating above physical ceiling
                        phys_max_caps.fetch_add(1, Ordering::Relaxed);
                        dlog_t = 0.0;
                    } else if t_old >= constants::T_HOT {
                        // Zone 2: Taper heating between T_HOT and T_PHYS_MAX
                        let scale = ((constants::T_PHYS_MAX - t_old) / (constants::T_PHYS_MAX - constants::T_HOT))
                            .clamp(0.0, 1.0);

                        hot_zone_entries.fetch_add(1, Ordering::Relaxed);
                        dlog_t *= scale;
                    }
                    // Zone 1 (t_old < T_HOT): Full compression heating (no change to dlog_t)
                }

                // Step 4: Apply exponential update
                let t_new = t_old * dlog_t.exp();

                // Hard panic on non-finite (catches numerical explosion early)
                if !t_new.is_finite() {
                    panic!(
                        "Non-finite temperature at ({}, {}): T_old={:.2} K, \
                         div_u={:.2} s⁻¹, dlog_t_raw={:.6}, dlog_t={:.6}",
                        x_u32, y_u32, t_old, div_u, dlog_t_raw, dlog_t
                    );
                }

                // Physical bounds (CMB minimum to physics ceiling)
                row[x].temperature = t_new.clamp(constants::T_CMB, constants::T_PHYS_MAX);
            }
        });

    // Log hot-spot brake activity (helps tune constants)
    let hot_entries = hot_zone_entries.load(Ordering::Relaxed);
    let phys_caps = phys_max_caps.load(Ordering::Relaxed);

    if hot_entries > 0 || phys_caps > 0 {
        bevy::log::debug!(
            hot_zone_entries = hot_entries,
            phys_max_caps = phys_caps,
            "Hot-spot brake activity this substep"
        );
    }
}

/// Compute velocity divergence at a given cell using central differences.
///
/// ∇·u = ∂u/∂x + ∂v/∂y
///     = (u_right - u_left)/(2dx) + (v_up - v_down)/(2dy)
///
/// Boundary conditions:
/// - At walls: use no-slip (u=0, v=0) for boundary velocity
/// - At domain edges: use no-slip (closed room assumption)
fn compute_divergence(
    cells: &[AtmosphereCell],
    width: u32,
    height: u32,
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    dx: f32,
) -> f32 {
    // X-direction velocity gradient: ∂u/∂x
    let u_right = if x < width - 1 && !collision_map.is_blocked(x + 1, y) {
        let idx = ((y * width) + (x + 1)) as usize;
        cells[idx].u
    } else {
        0.0 // Wall or boundary: no-slip
    };

    let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
        let idx = ((y * width) + (x - 1)) as usize;
        cells[idx].u
    } else {
        0.0 // Wall or boundary: no-slip
    };

    let du_dx = (u_right - u_left) / (2.0 * dx);

    // Y-direction velocity gradient: ∂v/∂y
    let v_up = if y < height - 1 && !collision_map.is_blocked(x, y + 1) {
        let idx = (((y + 1) * width) + x) as usize;
        cells[idx].v
    } else {
        0.0 // Wall or boundary: no-slip
    };

    let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
        let idx = (((y - 1) * width) + x) as usize;
        cells[idx].v
    } else {
        0.0 // Wall or boundary: no-slip
    };

    let dv_dy = (v_up - v_down) / (2.0 * dx);

    // Divergence: ∇·u = ∂u/∂x + ∂v/∂y
    du_dx + dv_dy
}

use crate::atmosphere::grid::AtmosphereCell;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expansion_cooling() {
        // Test that expanding gas (∇·u > 0) cools down
        let mut grid = AtmosphereGrid::new(3, 3, 1.0, 1.0);
        let collision_map = TileCollisionMap::test_map(3, 3);

        // Set up expanding velocity field: u increases to the right
        for y in 0..3 {
            for x in 0..3 {
                let idx = grid.index(x, y);
                grid.cells[idx].u = x as f32; // Velocity increases → ∂u/∂x > 0
                grid.cells[idx].v = 0.0;
                grid.cells[idx].temperature = 300.0;
            }
        }

        let temp_before = grid.cells[grid.index(1, 1)].temperature;

        compression_heating_step(&mut grid, &collision_map, 0.01);

        let temp_after = grid.cells[grid.index(1, 1)].temperature;

        // Expanding flow should cool the gas
        assert!(
            temp_after < temp_before,
            "Expansion (∇·u > 0) should cool gas, but {} >= {}",
            temp_after,
            temp_before
        );
    }

    #[test]
    fn test_compression_heating() {
        // Test that compressing gas (∇·u < 0) heats up
        let mut grid = AtmosphereGrid::new(3, 3, 1.0, 1.0);
        let collision_map = TileCollisionMap::test_map(3, 3);

        // Set up compressing velocity field: u decreases to the right
        for y in 0..3 {
            for x in 0..3 {
                let idx = grid.index(x, y);
                grid.cells[idx].u = -(x as f32); // Velocity decreases → ∂u/∂x < 0
                grid.cells[idx].v = 0.0;
                grid.cells[idx].temperature = 300.0;
            }
        }

        let temp_before = grid.cells[grid.index(1, 1)].temperature;

        compression_heating_step(&mut grid, &collision_map, 0.01);

        let temp_after = grid.cells[grid.index(1, 1)].temperature;

        // Compressing flow should heat the gas
        assert!(
            temp_after > temp_before,
            "Compression (∇·u < 0) should heat gas, but {} <= {}",
            temp_after,
            temp_before
        );
    }
}
