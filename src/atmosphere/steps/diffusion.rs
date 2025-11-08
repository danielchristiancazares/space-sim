use crate::atmosphere::constants;
use crate::atmosphere::grid::{AtmosphereCell, AtmosphereGrid};
use crate::tilemap::TileCollisionMap;

pub fn diffusion_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    let dx = atmosphere.tile_size_physical;

    const MAX_ALPHA: f32 = 0.2; // Safety factor below theoretical limit of 0.25

    let alpha_momentum = (constants::KINEMATIC_VISCOSITY * dt / (dx * dx)).min(MAX_ALPHA);
    let alpha_temperature = (constants::THERMAL_DIFFUSIVITY * dt / (dx * dx)).min(MAX_ALPHA);
    let d_coeff = (constants::GAS_DIFFUSIVITY * dt / (dx * dx)).min(MAX_ALPHA);

    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);

            let c = &atmosphere.cells_buffer[idx];
            let left = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32 - 1,
                y as i32,
                c,
            );
            let right = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32 + 1,
                y as i32,
                c,
            );
            let down = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32,
                y as i32 - 1,
                c,
            );
            let up = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32,
                y as i32 + 1,
                c,
            );

            let laplacian_u = left.u + right.u + down.u + up.u - 4.0 * c.u;
            let laplacian_v = left.v + right.v + down.v + up.v - 4.0 * c.v;

            atmosphere.cells[idx].u += alpha_momentum * laplacian_u;
            atmosphere.cells[idx].v += alpha_momentum * laplacian_v;

            let laplacian_t =
                left.temperature + right.temperature + down.temperature + up.temperature
                    - 4.0 * c.temperature;
            atmosphere.cells[idx].temperature += alpha_temperature * laplacian_t;

            let laplacian_o2 =
                left.rho_o2 + right.rho_o2 + down.rho_o2 + up.rho_o2 - 4.0 * c.rho_o2;
            let laplacian_n2 =
                left.rho_n2 + right.rho_n2 + down.rho_n2 + up.rho_n2 - 4.0 * c.rho_n2;
            let laplacian_co2 =
                left.rho_co2 + right.rho_co2 + down.rho_co2 + up.rho_co2 - 4.0 * c.rho_co2;

            atmosphere.cells[idx].rho_o2 += d_coeff * laplacian_o2;
            atmosphere.cells[idx].rho_n2 += d_coeff * laplacian_n2;
            atmosphere.cells[idx].rho_co2 += d_coeff * laplacian_co2;
        }
    }
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
}
