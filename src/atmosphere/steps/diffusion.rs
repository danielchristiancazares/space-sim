use crate::atmosphere::constants;
use crate::atmosphere::grid::{AtmosphereCell, AtmosphereGrid};
use crate::tilemap::TileCollisionMap;
use rayon::prelude::*;

/// Compute grid-scale artificial viscosity for shock regularization.
///
/// At coarse grid resolution (1m cells), molecular viscosity alone is insufficient
/// to regularize supersonic expansion shocks. Artificial viscosity represents
/// unresolved subgrid turbulence and prevents numerical oscillations.
///
/// Formula: ν_artificial = COEFF × Δx × |u_max|
///
/// Returns: Artificial kinematic viscosity [m²/s]
fn compute_artificial_viscosity(atmosphere: &AtmosphereGrid, dx: f32) -> f32 {
    let max_velocity = atmosphere
        .cells
        .iter()
        .map(|c| (c.u * c.u + c.v * c.v).sqrt())
        .fold(0.0, f32::max);

    constants::ARTIFICIAL_VISCOSITY * dx * max_velocity
}

pub fn diffusion_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    let dx = atmosphere.tile_size_physical;

    const MAX_ALPHA: f32 = 0.2; // Safety factor below theoretical limit of 0.25

    // Total viscosity: molecular + artificial (grid-scale regularization)
    let nu_molecular = constants::KINEMATIC_VISCOSITY; // 1.5e-5 m²/s (air at 20°C)
    let nu_artificial = compute_artificial_viscosity(atmosphere, dx);
    let nu_total = nu_molecular + nu_artificial;

    let alpha_momentum = (nu_total * dt / (dx * dx)).min(MAX_ALPHA);
    let alpha_temperature = (constants::THERMAL_DIFFUSIVITY * dt / (dx * dx)).min(MAX_ALPHA);
    let d_coeff = (constants::GAS_DIFFUSIVITY * dt / (dx * dx)).min(MAX_ALPHA);

    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    let width = atmosphere.width as usize;
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

                let c = &cells_buffer[idx];
                let left = get_or_boundary(
                    cells_buffer,
                    atmosphere.width,
                    atmosphere.height,
                    collision_map,
                    x_u32 as i32 - 1,
                    y_u32 as i32,
                    c,
                );
                let right = get_or_boundary(
                    cells_buffer,
                    atmosphere.width,
                    atmosphere.height,
                    collision_map,
                    x_u32 as i32 + 1,
                    y_u32 as i32,
                    c,
                );
                let down = get_or_boundary(
                    cells_buffer,
                    atmosphere.width,
                    atmosphere.height,
                    collision_map,
                    x_u32 as i32,
                    y_u32 as i32 - 1,
                    c,
                );
                let up = get_or_boundary(
                    cells_buffer,
                    atmosphere.width,
                    atmosphere.height,
                    collision_map,
                    x_u32 as i32,
                    y_u32 as i32 + 1,
                    c,
                );

                let laplacian_u = left.u + right.u + down.u + up.u - 4.0 * c.u;
                let laplacian_v = left.v + right.v + down.v + up.v - 4.0 * c.v;

                let cell = &mut row[x];
                cell.u += alpha_momentum * laplacian_u;
                cell.v += alpha_momentum * laplacian_v;

                let laplacian_t =
                    left.temperature + right.temperature + down.temperature + up.temperature
                        - 4.0 * c.temperature;
                cell.temperature += alpha_temperature * laplacian_t;

                let laplacian_o2 =
                    left.rho_o2 + right.rho_o2 + down.rho_o2 + up.rho_o2 - 4.0 * c.rho_o2;
                let laplacian_n2 =
                    left.rho_n2 + right.rho_n2 + down.rho_n2 + up.rho_n2 - 4.0 * c.rho_n2;
                let laplacian_co2 =
                    left.rho_co2 + right.rho_co2 + down.rho_co2 + up.rho_co2 - 4.0 * c.rho_co2;

                cell.rho_o2 += d_coeff * laplacian_o2;
                cell.rho_n2 += d_coeff * laplacian_n2;
                cell.rho_co2 += d_coeff * laplacian_co2;
            }
        });
}

fn get_or_boundary(
    cells: &[AtmosphereCell],
    width: u32,
    height: u32,
    collision_map: &TileCollisionMap,
    x: i32,
    y: i32,
    center: &AtmosphereCell,
) -> AtmosphereCell {
    if x >= 0 && y >= 0 && (x as u32) < width && (y as u32) < height {
        let idx = (y as u32 * width + x as u32) as usize;
        if !collision_map.is_blocked(x as u32, y as u32) {
            cells[idx].clone()
        } else {
            AtmosphereCell {
                rho_o2: center.rho_o2,
                rho_n2: center.rho_n2,
                rho_co2: center.rho_co2,
                u: 0.0,
                v: 0.0,
                temperature: center.temperature,
                pressure: center.pressure,
            }
        }
    } else {
        // Treat outside-the-grid as solid walls with no-slip boundary condition.
        // Zero-gradient for scalars (mass, temperature), no-slip (u=0, v=0) for velocity.
        // This prevents mass diffusion while enforcing physical boundary conditions.
        AtmosphereCell {
            rho_o2: center.rho_o2,
            rho_n2: center.rho_n2,
            rho_co2: center.rho_co2,
            u: 0.0, // No-slip: velocity must be zero at walls
            v: 0.0, // No-slip: velocity must be zero at walls
            temperature: center.temperature,
            pressure: center.pressure,
        }
    }
}
