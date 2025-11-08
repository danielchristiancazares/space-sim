use crate::atmosphere::grid::{AtmosphereCell, AtmosphereGrid};
use crate::tilemap::TileCollisionMap;

pub fn advection_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    // Copy current state to buffer (original state φⁿ)
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    // Step 1: Forward advection (predictor)
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells_buffer[idx];

            let x_back = x as f32 - (cell.u * dt) / atmosphere.tile_size_physical;
            let y_back = y as f32 - (cell.v * dt) / atmosphere.tile_size_physical;

            atmosphere.forward_state[idx] = interpolate_cell(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                x_back,
                y_back,
            );
        }
    }

    // Step 2: Backward advection (corrector)
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.forward_state[idx];

            let x_forward = x as f32 + (cell.u * dt) / atmosphere.tile_size_physical;
            let y_forward = y as f32 + (cell.v * dt) / atmosphere.tile_size_physical;

            atmosphere.backward_state[idx] = interpolate_cell(
                &atmosphere.forward_state,
                atmosphere.width,
                atmosphere.height,
                x_forward,
                y_forward,
            );
        }
    }

    // Step 3: MacCormack correction
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let original = &atmosphere.cells_buffer[idx];
            let forward = &atmosphere.forward_state[idx];
            let backward = &atmosphere.backward_state[idx];

            let corrected_rho_o2 = forward.rho_o2 + 0.5 * (original.rho_o2 - backward.rho_o2);
            let corrected_rho_n2 = forward.rho_n2 + 0.5 * (original.rho_n2 - backward.rho_n2);
            let corrected_rho_co2 = forward.rho_co2 + 0.5 * (original.rho_co2 - backward.rho_co2);
            let corrected_u = forward.u + 0.5 * (original.u - backward.u);
            let corrected_v = forward.v + 0.5 * (original.v - backward.v);
            let corrected_temp =
                forward.temperature + 0.5 * (original.temperature - backward.temperature);

            atmosphere.cells[idx].rho_o2 = corrected_rho_o2.max(0.0);
            atmosphere.cells[idx].rho_n2 = corrected_rho_n2.max(0.0);
            atmosphere.cells[idx].rho_co2 = corrected_rho_co2.max(0.0);
            atmosphere.cells[idx].u = corrected_u;
            atmosphere.cells[idx].v = corrected_v;
            atmosphere.cells[idx].temperature = corrected_temp.max(0.0);
            // Pressure updated globally between steps
        }
    }
}

fn interpolate_cell(
    cells: &[AtmosphereCell],
    width: u32,
    height: u32,
    x: f32,
    y: f32,
) -> AtmosphereCell {
    let x0 = x.floor().max(0.0).min((width - 1) as f32) as u32;
    let y0 = y.floor().max(0.0).min((height - 1) as f32) as u32;
    let x1 = (x0 + 1).min(width - 1);
    let y1 = (y0 + 1).min(height - 1);

    let fx = (x - x0 as f32).clamp(0.0, 1.0);
    let fy = (y - y0 as f32).clamp(0.0, 1.0);

    let idx00 = (y0 * width + x0) as usize;
    let idx10 = (y0 * width + x1) as usize;
    let idx01 = (y1 * width + x0) as usize;
    let idx11 = (y1 * width + x1) as usize;

    let c00 = &cells[idx00];
    let c10 = &cells[idx10];
    let c01 = &cells[idx01];
    let c11 = &cells[idx11];

    AtmosphereCell {
        rho_o2: bilerp(c00.rho_o2, c10.rho_o2, c01.rho_o2, c11.rho_o2, fx, fy),
        rho_n2: bilerp(c00.rho_n2, c10.rho_n2, c01.rho_n2, c11.rho_n2, fx, fy),
        rho_co2: bilerp(c00.rho_co2, c10.rho_co2, c01.rho_co2, c11.rho_co2, fx, fy),
        u: bilerp(c00.u, c10.u, c01.u, c11.u, fx, fy),
        v: bilerp(c00.v, c10.v, c01.v, c11.v, fx, fy),
        temperature: bilerp(
            c00.temperature,
            c10.temperature,
            c01.temperature,
            c11.temperature,
            fx,
            fy,
        ),
        pressure: 0.0,
    }
}

#[inline]
fn bilerp(v00: f32, v10: f32, v01: f32, v11: f32, fx: f32, fy: f32) -> f32 {
    let v0 = v00 * (1.0 - fx) + v10 * fx;
    let v1 = v01 * (1.0 - fx) + v11 * fx;
    v0 * (1.0 - fy) + v1 * fy
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tilemap::TileCollisionMap;

    #[test]
    fn test_maccormack_differs_from_semi_lagrangian() {
        let mut mac_grid = AtmosphereGrid::new(3, 1, 1.0, 1.0);
        let mut semi_grid = AtmosphereGrid::new(3, 1, 1.0, 1.0);

        // Create a sharp density gradient that will advect to the right.
        let setup = [
            AtmosphereCell {
                rho_o2: 1.0,
                rho_n2: 0.0,
                rho_co2: 0.0,
                u: 1.0,
                v: 0.0,
                temperature: 300.0,
                pressure: 0.0,
            },
            AtmosphereCell {
                rho_o2: 0.25,
                rho_n2: 0.0,
                rho_co2: 0.0,
                u: 1.0,
                v: 0.0,
                temperature: 300.0,
                pressure: 0.0,
            },
            AtmosphereCell {
                rho_o2: 0.0,
                rho_n2: 0.0,
                rho_co2: 0.0,
                u: 1.0,
                v: 0.0,
                temperature: 300.0,
                pressure: 0.0,
            },
        ];

        mac_grid.cells.clone_from_slice(&setup);
        semi_grid.cells.clone_from_slice(&setup);

        let collision_map = TileCollisionMap::test_map(3, 1);
        let dt = 0.4;

        advection_step(&mut mac_grid, &collision_map, dt);
        semi_lagrangian_step(&mut semi_grid, &collision_map, dt);

        let mac_density = mac_grid.cells[1].rho_o2;
        let semi_density = semi_grid.cells[1].rho_o2;

        assert!(
            (mac_density - semi_density).abs() > 1e-4,
            "MacCormack path should differ from semi-Lagrangian reference"
        );
    }

    #[test]
    fn test_bilinear_interpolation() {
        let v00 = 0.0;
        let v10 = 1.0;
        let v01 = 0.0;
        let v11 = 1.0;

        let center = bilerp(v00, v10, v01, v11, 0.5, 0.5);
        assert!((center - 0.5).abs() < 1e-6);

        let corner_00 = bilerp(v00, v10, v01, v11, 0.0, 0.0);
        assert!((corner_00 - v00).abs() < 1e-6);

        let corner_11 = bilerp(v00, v10, v01, v11, 1.0, 1.0);
        assert!((corner_11 - v11).abs() < 1e-6);
    }
}

/// Minimal semi-Lagrangian integrator used to compare against the MacCormack path.
#[cfg(test)]
fn semi_lagrangian_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    atmosphere
        .cells_buffer
        .clone_from_slice(&atmosphere.cells);

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells_buffer[idx];

            let x_back = x as f32 - (cell.u * dt) / atmosphere.tile_size_physical;
            let y_back = y as f32 - (cell.v * dt) / atmosphere.tile_size_physical;

            atmosphere.cells[idx] = interpolate_cell(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                x_back,
                y_back,
            );
        }
    }
}
