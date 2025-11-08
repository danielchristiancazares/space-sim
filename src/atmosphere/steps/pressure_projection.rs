use crate::atmosphere::constants;
use crate::atmosphere::grid::AtmosphereGrid;
use crate::tilemap::TileCollisionMap;

pub fn pressure_correction_step(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;

    for v in atmosphere.pressure_correction.iter_mut() {
        *v = 0.0;
    }
    for v in atmosphere.pressure_correction_temp.iter_mut() {
        *v = 0.0;
    }
    for v in atmosphere.divergence.iter_mut() {
        *v = 0.0;
    }

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

            atmosphere.divergence[idx] = du_dx + dv_dy;
        }
    }

    for _ in 0..constants::POISSON_MAX_ITERATIONS {
        let mut max_change: f32 = 0.0;

        for y in 0..atmosphere.height {
            for x in 0..atmosphere.width {
                if collision_map.is_blocked(x, y) {
                    continue;
                }

                let idx = atmosphere.index(x, y);
                let rho = atmosphere.cells[idx].total_density().max(1e-6);

                let p_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                    atmosphere.pressure_correction[atmosphere.index(x - 1, y)]
                } else {
                    atmosphere.pressure_correction[idx]
                };

                let p_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                    atmosphere.pressure_correction[atmosphere.index(x + 1, y)]
                } else {
                    atmosphere.pressure_correction[idx]
                };

                let p_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                    atmosphere.pressure_correction[atmosphere.index(x, y - 1)]
                } else {
                    atmosphere.pressure_correction[idx]
                };

                let p_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                    atmosphere.pressure_correction[atmosphere.index(x, y + 1)]
                } else {
                    atmosphere.pressure_correction[idx]
                };

                let rhs = (rho / dt) * atmosphere.divergence[idx];
                let p_new = (p_left + p_right + p_down + p_up - dx * dx * rhs) * 0.25;

                atmosphere.pressure_correction_temp[idx] = p_new;
                max_change = max_change.max((p_new - atmosphere.pressure_correction[idx]).abs());
            }
        }

        atmosphere
            .pressure_correction
            .copy_from_slice(&atmosphere.pressure_correction_temp);

        if max_change < constants::POISSON_TOLERANCE {
            break;
        }
    }

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let rho = atmosphere.cells[idx].total_density().max(1e-6);

            let p_c = atmosphere.pressure_correction[idx];

            let p_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                atmosphere.pressure_correction[atmosphere.index(x - 1, y)]
            } else {
                p_c
            };

            let p_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                atmosphere.pressure_correction[atmosphere.index(x + 1, y)]
            } else {
                p_c
            };

            let p_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                atmosphere.pressure_correction[atmosphere.index(x, y - 1)]
            } else {
                p_c
            };

            let p_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                atmosphere.pressure_correction[atmosphere.index(x, y + 1)]
            } else {
                p_c
            };

            let dp_dx = (p_right - p_left) / (2.0 * dx);
            let dp_dy = (p_up - p_down) / (2.0 * dx);

            atmosphere.cells[idx].u -= (dt / rho) * dp_dx;
            atmosphere.cells[idx].v -= (dt / rho) * dp_dy;
        }
    }
}
