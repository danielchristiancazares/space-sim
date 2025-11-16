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

    // Physical bounds for temperature
    const T_MIN: f32 = 2.7; // K (cosmic microwave background)
    const T_MAX: f32 = 10_000.0; // K (safety cap for numerical stability)

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

                let t_old = cell.temperature;
                let gamma_minus_1 = constants::GAMMA - 1.0;

                // ============================================================
                // Temperature compression: ∂T/∂t = -(γ-1)T∇·u (energy equation)
                // ============================================================
                // EXACT EXPONENTIAL UPDATE:
                // Treating ∇·u as constant over substep, the exact solution is:
                //   T(t+dt) = T(t) × exp[-(γ-1) × ∇·u × dt]
                //
                // This is mathematically precise and prevents compounding errors.
                let exponent = -gamma_minus_1 * div_u * dt;

                // Guard exponent to catch divergence × dt explosions
                debug_assert!(
                    exponent.abs() < 10.0,
                    "Extreme compression exponent at ({}, {}): exp({:.3}) from \
                     div_u={:.2} s⁻¹, dt={:.4} s (CFL violation likely)",
                    x_u32,
                    y_u32,
                    exponent,
                    div_u,
                    dt
                );

                let factor = exponent.exp();
                let t_new = t_old * factor;

                // Hard panic on non-finite (catches numerical explosion early)
                if !t_new.is_finite() {
                    panic!(
                        "Non-finite temperature at ({}, {}): T_old={:.2} K, \
                         div_u={:.2} s⁻¹, exponent={:.3}, factor={}",
                        x_u32, y_u32, t_old, div_u, exponent, factor
                    );
                }

                // Physical bounds (CMB minimum to numerical safety cap)
                row[x].temperature = t_new.clamp(T_MIN, T_MAX);
            }
        });
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
