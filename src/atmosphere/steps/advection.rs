use crate::atmosphere::grid::{AtmosphereCell, AtmosphereGrid};
use crate::tilemap::TileCollisionMap;
use rayon::prelude::*;

pub fn advection_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    // Copy current state to buffer (original state φⁿ)
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    let _track_negatives = cfg!(debug_assertions); // No longer used (was for MacCormack debugging)
    let _tile_volume = atmosphere.tile_size_physical
        * atmosphere.tile_size_physical
        * crate::atmosphere::constants::ROOM_HEIGHT;

    // STEP 1: Flux-conservative density advection
    // Directly discretizes: ∂ρ/∂t = -∇·(ρu)
    advect_density_flux_conservative(atmosphere, collision_map, dt);

    // STEP 2: Rusanov flux advection for velocity and temperature
    // Upwind monotone scheme prevents oscillations at shocks
    advect_velocity_temperature_rusanov(atmosphere, collision_map, dt);
}

/// Flux-conservative density advection using 1st order upwind scheme.
///
/// Directly discretizes: ∂ρ/∂t = -∇·(ρu)
///
/// Guarantees mass conservation by explicitly tracking flux across cell faces:
/// - flux out of cell A = flux into cell B
/// - No interpolation artifacts
/// - Conservative by construction
fn advect_density_flux_conservative(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;
    let dy = dx; // Square cells
    let cell_area = dx * dy;

    let cells_buffer = &atmosphere.cells_buffer;
    let width = atmosphere.width;
    let height = atmosphere.height;

    for y in 0..height {
        for x in 0..width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y) as usize;

            // Compute signed fluxes at each face
            let (flux_e_o2, flux_e_n2, flux_e_co2) =
                compute_flux_east(cells_buffer, collision_map, x, y, width, height, dy);
            let (flux_w_o2, flux_w_n2, flux_w_co2) =
                compute_flux_west(cells_buffer, collision_map, x, y, width, height, dy);
            let (flux_n_o2, flux_n_n2, flux_n_co2) =
                compute_flux_north(cells_buffer, collision_map, x, y, width, height, dx);
            let (flux_s_o2, flux_s_n2, flux_s_co2) =
                compute_flux_south(cells_buffer, collision_map, x, y, width, height, dx);

            // Net outward flux
            let net_flux_out_o2 = flux_e_o2 - flux_w_o2 + flux_n_o2 - flux_s_o2;
            let net_flux_out_n2 = flux_e_n2 - flux_w_n2 + flux_n_n2 - flux_s_n2;
            let net_flux_out_co2 = flux_e_co2 - flux_w_co2 + flux_n_co2 - flux_s_co2;

            // Update densities: ∂ρ/∂t = -net_flux_out / area
            atmosphere.cells[idx].rho_o2 -= (net_flux_out_o2 / cell_area) * dt;
            atmosphere.cells[idx].rho_n2 -= (net_flux_out_n2 / cell_area) * dt;
            atmosphere.cells[idx].rho_co2 -= (net_flux_out_co2 / cell_area) * dt;

            // Clamp to physical minimum
            atmosphere.cells[idx].rho_o2 = atmosphere.cells[idx].rho_o2.max(1e-10);
            atmosphere.cells[idx].rho_n2 = atmosphere.cells[idx].rho_n2.max(1e-10);
            atmosphere.cells[idx].rho_co2 = atmosphere.cells[idx].rho_co2.max(1e-10);
        }
    }
}

/// Compute flux at east face (x+1/2) for all species.
/// Returns: (flux_o2, flux_n2, flux_co2) in [kg/s]
fn compute_flux_east(
    cells: &[AtmosphereCell],
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    width: u32,
    _height: u32,
    dy: f32,
) -> (f32, f32, f32) {
    if x >= width - 1 || collision_map.is_blocked(x + 1, y) {
        return (0.0, 0.0, 0.0); // No flux at boundary/wall
    }

    let idx_l = (y * width + x) as usize;
    let idx_r = (y * width + (x + 1)) as usize;

    // Face velocity (average)
    let u_face = 0.5 * (cells[idx_l].u + cells[idx_r].u);

    // Upwind density selection
    let (rho_o2, rho_n2, rho_co2) = if u_face > 0.0 {
        (
            cells[idx_l].rho_o2,
            cells[idx_l].rho_n2,
            cells[idx_l].rho_co2,
        ) // Flow →
    } else {
        (
            cells[idx_r].rho_o2,
            cells[idx_r].rho_n2,
            cells[idx_r].rho_co2,
        ) // Flow ←
    };

    // Flux = ρ × u × area (kg/s)
    let flux_o2 = rho_o2 * u_face * dy;
    let flux_n2 = rho_n2 * u_face * dy;
    let flux_co2 = rho_co2 * u_face * dy;

    (flux_o2, flux_n2, flux_co2)
}

/// Compute flux at west face (x-1/2) for all species.
fn compute_flux_west(
    cells: &[AtmosphereCell],
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    width: u32,
    _height: u32,
    dy: f32,
) -> (f32, f32, f32) {
    if x == 0 || collision_map.is_blocked(x - 1, y) {
        return (0.0, 0.0, 0.0);
    }

    let idx_l = (y * width + (x - 1)) as usize;
    let idx_r = (y * width + x) as usize;

    let u_face = 0.5 * (cells[idx_l].u + cells[idx_r].u);

    let (rho_o2, rho_n2, rho_co2) = if u_face > 0.0 {
        (
            cells[idx_l].rho_o2,
            cells[idx_l].rho_n2,
            cells[idx_l].rho_co2,
        )
    } else {
        (
            cells[idx_r].rho_o2,
            cells[idx_r].rho_n2,
            cells[idx_r].rho_co2,
        )
    };

    (
        rho_o2 * u_face * dy,
        rho_n2 * u_face * dy,
        rho_co2 * u_face * dy,
    )
}

/// Compute flux at north face (y+1/2) for all species.
fn compute_flux_north(
    cells: &[AtmosphereCell],
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    width: u32,
    height: u32,
    dx: f32,
) -> (f32, f32, f32) {
    if y >= height - 1 || collision_map.is_blocked(x, y + 1) {
        return (0.0, 0.0, 0.0);
    }

    let idx_s = (y * width + x) as usize;
    let idx_n = ((y + 1) * width + x) as usize;

    let v_face = 0.5 * (cells[idx_s].v + cells[idx_n].v);

    let (rho_o2, rho_n2, rho_co2) = if v_face > 0.0 {
        (
            cells[idx_s].rho_o2,
            cells[idx_s].rho_n2,
            cells[idx_s].rho_co2,
        ) // Flow ↑
    } else {
        (
            cells[idx_n].rho_o2,
            cells[idx_n].rho_n2,
            cells[idx_n].rho_co2,
        ) // Flow ↓
    };

    (
        rho_o2 * v_face * dx,
        rho_n2 * v_face * dx,
        rho_co2 * v_face * dx,
    )
}

/// Compute flux at south face (y-1/2) for all species.
fn compute_flux_south(
    cells: &[AtmosphereCell],
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    width: u32,
    _height: u32,
    dx: f32,
) -> (f32, f32, f32) {
    if y == 0 || collision_map.is_blocked(x, y - 1) {
        return (0.0, 0.0, 0.0);
    }

    let idx_s = ((y - 1) * width + x) as usize;
    let idx_n = (y * width + x) as usize;

    let v_face = 0.5 * (cells[idx_s].v + cells[idx_n].v);

    let (rho_o2, rho_n2, rho_co2) = if v_face > 0.0 {
        (
            cells[idx_s].rho_o2,
            cells[idx_s].rho_n2,
            cells[idx_s].rho_co2,
        )
    } else {
        (
            cells[idx_n].rho_o2,
            cells[idx_n].rho_n2,
            cells[idx_n].rho_co2,
        )
    };

    (
        rho_o2 * v_face * dx,
        rho_n2 * v_face * dx,
        rho_co2 * v_face * dx,
    )
}

/// Rusanov (local Lax-Friedrichs) flux for velocity and temperature advection.
///
/// Implements monotone upwind scheme with numerical dissipation:
///   flux = 0.5 * [F(q_L) + F(q_R) - |a| * (q_R - q_L)]
///
/// Where:
/// - F(q) = u * q (physical flux)  
/// - a = |u_face| (local wave speed for convective transport)
/// - Upwind dissipation term ensures monotonicity (no oscillations at shocks)
///
/// This is 1st-order accurate but unconditionally monotone, ideal for
/// capturing shocks and contact discontinuities without spurious oscillations.
fn advect_velocity_temperature_rusanov(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;
    let width = atmosphere.width as usize;

    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);
    let cells = &atmosphere.cells_buffer;

    // Parallel update over rows
    atmosphere
        .cells
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

                // East face (x+1/2)
                let (flux_u_e, flux_v_e, flux_t_e) = if x_u32 < atmosphere.width - 1
                    && !collision_map.is_blocked(x_u32 + 1, y_u32)
                {
                    let idx_r = y * width + (x + 1);

                    let q_l = cells[idx].u;
                    let q_r = cells[idx_r].u;
                    let u_face = 0.5 * (cells[idx].u + cells[idx_r].u);
                    let f_l = u_face * q_l;
                    let f_r = u_face * q_r;
                    let a = u_face.abs();
                    let flux_u = 0.5 * (f_l + f_r - a * (q_r - q_l));

                    let q_l_v = cells[idx].v;
                    let q_r_v = cells[idx_r].v;
                    let f_l_v = u_face * q_l_v;
                    let f_r_v = u_face * q_r_v;
                    let flux_v = 0.5 * (f_l_v + f_r_v - a * (q_r_v - q_l_v));

                    let q_l_t = cells[idx].temperature;
                    let q_r_t = cells[idx_r].temperature;
                    let f_l_t = u_face * q_l_t;
                    let f_r_t = u_face * q_r_t;
                    let flux_t = 0.5 * (f_l_t + f_r_t - a * (q_r_t - q_l_t));

                    (flux_u, flux_v, flux_t)
                } else {
                    (0.0, 0.0, 0.0)
                };

                // West face (x-1/2)
                let (flux_u_w, flux_v_w, flux_t_w) =
                    if x_u32 > 0 && !collision_map.is_blocked(x_u32 - 1, y_u32) {
                        let idx_l = y * width + (x - 1);

                        let q_l = cells[idx_l].u;
                        let q_r = cells[idx].u;
                        let u_face = 0.5 * (cells[idx_l].u + cells[idx].u);
                        let f_l = u_face * q_l;
                        let f_r = u_face * q_r;
                        let a = u_face.abs();
                        let flux_u = 0.5 * (f_l + f_r - a * (q_r - q_l));

                        let q_l_v = cells[idx_l].v;
                        let q_r_v = cells[idx].v;
                        let f_l_v = u_face * q_l_v;
                        let f_r_v = u_face * q_r_v;
                        let flux_v = 0.5 * (f_l_v + f_r_v - a * (q_r_v - q_l_v));

                        let q_l_t = cells[idx_l].temperature;
                        let q_r_t = cells[idx].temperature;
                        let f_l_t = u_face * q_l_t;
                        let f_r_t = u_face * q_r_t;
                        let flux_t = 0.5 * (f_l_t + f_r_t - a * (q_r_t - q_l_t));

                        (flux_u, flux_v, flux_t)
                    } else {
                        (0.0, 0.0, 0.0)
                    };

                // North face (y+1/2)
                let (flux_u_n, flux_v_n, flux_t_n) = if y_u32 < atmosphere.height - 1
                    && !collision_map.is_blocked(x_u32, y_u32 + 1)
                {
                    let idx_n = (y + 1) * width + x;

                    let q_s = cells[idx].u;
                    let q_n = cells[idx_n].u;
                    let v_face = 0.5 * (cells[idx].v + cells[idx_n].v);
                    let f_s = v_face * q_s;
                    let f_n = v_face * q_n;
                    let a = v_face.abs();
                    let flux_u = 0.5 * (f_s + f_n - a * (q_n - q_s));

                    let q_s_v = cells[idx].v;
                    let q_n_v = cells[idx_n].v;
                    let f_s_v = v_face * q_s_v;
                    let f_n_v = v_face * q_n_v;
                    let flux_v = 0.5 * (f_s_v + f_n_v - a * (q_n_v - q_s_v));

                    let q_s_t = cells[idx].temperature;
                    let q_n_t = cells[idx_n].temperature;
                    let f_s_t = v_face * q_s_t;
                    let f_n_t = v_face * q_n_t;
                    let flux_t = 0.5 * (f_s_t + f_n_t - a * (q_n_t - q_s_t));

                    (flux_u, flux_v, flux_t)
                } else {
                    (0.0, 0.0, 0.0)
                };

                // South face (y-1/2)
                let (flux_u_s, flux_v_s, flux_t_s) =
                    if y_u32 > 0 && !collision_map.is_blocked(x_u32, y_u32 - 1) {
                        let idx_s = (y - 1) * width + x;

                        let q_s = cells[idx_s].u;
                        let q_n = cells[idx].u;
                        let v_face = 0.5 * (cells[idx_s].v + cells[idx].v);
                        let f_s = v_face * q_s;
                        let f_n = v_face * q_n;
                        let a = v_face.abs();
                        let flux_u = 0.5 * (f_s + f_n - a * (q_n - q_s));

                        let q_s_v = cells[idx_s].v;
                        let q_n_v = cells[idx].v;
                        let f_s_v = v_face * q_s_v;
                        let f_n_v = v_face * q_n_v;
                        let flux_v = 0.5 * (f_s_v + f_n_v - a * (q_n_v - q_s_v));

                        let q_s_t = cells[idx_s].temperature;
                        let q_n_t = cells[idx].temperature;
                        let f_s_t = v_face * q_s_t;
                        let f_n_t = v_face * q_n_t;
                        let flux_t = 0.5 * (f_s_t + f_n_t - a * (q_n_t - q_s_t));

                        (flux_u, flux_v, flux_t)
                    } else {
                        (0.0, 0.0, 0.0)
                    };

                // Update: dq/dt = -(flux_E - flux_W + flux_N - flux_S) / dx
                row[x].u -= (flux_u_e - flux_u_w + flux_u_n - flux_u_s) * dt / dx;
                row[x].v -= (flux_v_e - flux_v_w + flux_v_n - flux_v_s) * dt / dx;
                row[x].temperature -= (flux_t_e - flux_t_w + flux_t_n - flux_t_s) * dt / dx;

                // Clamp temperature to CMB minimum
                use crate::atmosphere::constants;
                row[x].temperature = row[x].temperature.max(constants::T_CMB);
            }
        });

    // Runtime guard for non-finite values (runs in --release)
    for (idx, cell) in atmosphere.cells.iter().enumerate() {
        if !cell.u.is_finite() || !cell.v.is_finite() || !cell.temperature.is_finite() {
            let x = idx % width;
            let y = idx / width;
            panic!(
                "Non-finite values after Rusanov at ({}, {}): u={}, v={}, T={}",
                x, y, cell.u, cell.v, cell.temperature
            );
        }
    }
}
