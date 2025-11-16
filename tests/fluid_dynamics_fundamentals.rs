/// Integration tests validating fundamental fluid dynamics behavior.
/// These tests check if the simulation produces physically plausible results
/// ("sanity checks") rather than high-precision accuracy.
use bevy::log::{Level, LogPlugin};
use bevy::prelude::*;
use bevy::time::TimeUpdateStrategy;
use fluid_dynamics_research_sim::{
    atmosphere::{constants, update_grid_pressures, AtmosphereCell, AtmosphereGrid},
    tilemap::TileCollisionMap,
};
use std::time::Duration;

const GRID_SIZE: u32 = 16; // Larger grid for fluid dynamics tests
const TILE_SIZE: f32 = 1.0;
const DT: f32 = 0.01; // Small timestep for stability
const MOMENTUM_ERROR_LIMIT: f32 = 0.05; // 5% relative drift tolerance
const MAX_DIVERGENCE_LIMIT: f32 = 0.5; // Allow small numerical divergence

/// Test that diffusion smooths out sharp gradients over time.
/// Creates a temperature spike and verifies it diffuses.
#[test]
fn diffusion_smooths_gradients() {
    let mut grid = AtmosphereGrid::new(GRID_SIZE, GRID_SIZE, TILE_SIZE, TILE_SIZE);
    let collision_map = TileCollisionMap::empty(GRID_SIZE, GRID_SIZE, Vec2::ONE, Vec2::ZERO);

    // Initialize: uniform density, temperature spike in center
    for y in 0..GRID_SIZE {
        for x in 0..GRID_SIZE {
            let idx = grid.index(x, y);
            let mut cell = AtmosphereCell::earth_atmosphere();
            cell.u = 0.0;
            cell.v = 0.0;

            if x == GRID_SIZE / 2 && y == GRID_SIZE / 2 {
                cell.temperature = 400.0; // Hot spot
            } else {
                cell.temperature = 293.0; // Room temp
            }

            grid.cells[idx] = cell;
        }
    }

    let center_idx = grid.index(GRID_SIZE / 2, GRID_SIZE / 2);
    let neighbor_idx = grid.index(GRID_SIZE / 2 + 1, GRID_SIZE / 2);

    let initial_center_temp = grid.cells[center_idx].temperature;
    let initial_neighbor_temp = grid.cells[neighbor_idx].temperature;
    let initial_gradient = initial_center_temp - initial_neighbor_temp;

    // Run diffusion steps
    for _ in 0..200 {
        fluid_dynamics_research_sim::atmosphere::steps::diffusion::diffusion_step(
            &mut grid,
            &collision_map,
            DT,
        );
    }

    let final_center_temp = grid.cells[center_idx].temperature;
    let final_neighbor_temp = grid.cells[neighbor_idx].temperature;
    let final_gradient = final_center_temp - final_neighbor_temp;

    // Verify: gradient should be smaller (diffusion smooths)
    assert!(
        final_gradient < initial_gradient,
        "Temperature gradient should reduce: {} K -> {} K",
        initial_gradient,
        final_gradient
    );

    // Verify: neighbor temperature should increase (heat spread)
    assert!(
        final_neighbor_temp > initial_neighbor_temp,
        "Neighbor should warm up: {} K -> {} K",
        initial_neighbor_temp,
        final_neighbor_temp
    );
}

/// Test that advection transports quantities along the velocity field.
/// Creates a uniform velocity field and a CO2 blob, checks that blob moves.
#[test]
fn advection_transports_tracers() {
    let mut grid = AtmosphereGrid::new(GRID_SIZE, GRID_SIZE, TILE_SIZE, TILE_SIZE);
    let collision_map = TileCollisionMap::empty(GRID_SIZE, GRID_SIZE, Vec2::ONE, Vec2::ZERO);

    // Initialize: uniform rightward velocity, CO2 blob on left
    for y in 0..GRID_SIZE {
        for x in 0..GRID_SIZE {
            let idx = grid.index(x, y);
            let mut cell = AtmosphereCell::earth_atmosphere();
            cell.u = 5.0; // Rightward velocity
            cell.v = 0.0;

            // CO2 tracer blob in left region
            if x < GRID_SIZE / 4 && y >= GRID_SIZE / 3 && y <= 2 * GRID_SIZE / 3 {
                cell.rho_co2 *= 10.0; // 10x normal CO2
            }

            grid.cells[idx] = cell;
        }
    }

    // Compute initial CO2 center of mass
    let initial_com = compute_co2_center_of_mass(&grid);

    // Run advection steps (should move CO2 blob rightward)
    for _ in 0..50 {
        fluid_dynamics_research_sim::atmosphere::steps::advection::advection_step(
            &mut grid,
            &collision_map,
            DT,
        );
    }

    let final_com = compute_co2_center_of_mass(&grid);

    // Verify: CO2 center of mass should move in direction of velocity
    let displacement_x = final_com.0 - initial_com.0;
    assert!(
        displacement_x > 0.5,
        "CO2 blob should move rightward: x_com {} -> {} (displacement: {})",
        initial_com.0,
        final_com.0,
        displacement_x
    );
}

/// Test that conservation laws hold over long simulation runs.
/// Runs simulation for many steps and verifies mass/momentum don't drift.
#[test]
fn conservation_long_term() {
    let mut app = App::new();
    app.add_plugins(MinimalPlugins);
    app.add_plugins(LogPlugin {
        level: Level::WARN, // Only warnings/errors to keep test output clean
        filter: "fluid_dynamics_research_sim=warn".into(),
        ..default()
    });

    let collision_map = TileCollisionMap::empty(GRID_SIZE, GRID_SIZE, Vec2::ONE, Vec2::ZERO);
    app.insert_resource(collision_map);

    let mut grid = AtmosphereGrid::new(GRID_SIZE, GRID_SIZE, TILE_SIZE, TILE_SIZE);

    // Initialize with Earth atmosphere everywhere
    for cell in grid.cells.iter_mut() {
        *cell = AtmosphereCell::earth_atmosphere();
    }

    app.insert_resource(grid);
    app.world_mut()
        .insert_resource(TimeUpdateStrategy::ManualDuration(Duration::from_secs_f32(
            DT,
        )));

    app.add_systems(Update, simulate_step.chain());

    // Compute initial conserved quantities
    let (initial_mass, initial_momentum) = {
        let world = app.world();
        let grid = world.resource::<AtmosphereGrid>();
        (compute_total_mass(grid), compute_total_momentum(grid))
    };

    // Run for many steps (10 simulated seconds)
    const STEPS: usize = 1000;
    for _ in 0..STEPS {
        app.update();
    }

    // Compute final conserved quantities
    let world = app.world();
    let grid = world.resource::<AtmosphereGrid>();
    let collision_map = world.resource::<TileCollisionMap>();
    let final_mass = compute_total_mass(grid);
    let final_momentum = compute_total_momentum(grid);

    // Verify: mass conservation (no sources/sinks, so should be exact)
    let mass_error = ((final_mass - initial_mass) / initial_mass).abs();
    assert!(
        mass_error < 0.001,
        "Mass should be conserved: initial={:.6} kg, final={:.6} kg, error={:.3}%",
        initial_mass,
        final_mass,
        mass_error * 100.0
    );

    // Verify: momentum conservation (no external forces)
    // Note: momentum can change due to pressure forces, but shouldn't drift systematically
    let momentum_x_error =
        ((final_momentum.0 - initial_momentum.0) / initial_momentum.0.abs().max(1.0)).abs();
    let momentum_y_error =
        ((final_momentum.1 - initial_momentum.1) / initial_momentum.1.abs().max(1.0)).abs();

    assert!(
        momentum_x_error < MOMENTUM_ERROR_LIMIT,
        "X-momentum drifted {:.3}% (allowed {:.3}%): initial={:.6}, final={:.6}",
        momentum_x_error * 100.0,
        MOMENTUM_ERROR_LIMIT * 100.0,
        initial_momentum.0,
        final_momentum.0
    );
    assert!(
        momentum_y_error < MOMENTUM_ERROR_LIMIT,
        "Y-momentum drifted {:.3}% (allowed {:.3}%): initial={:.6}, final={:.6}",
        momentum_y_error * 100.0,
        MOMENTUM_ERROR_LIMIT * 100.0,
        initial_momentum.1,
        final_momentum.1
    );

    let max_divergence = compute_max_divergence(grid, &collision_map);
    assert!(
        max_divergence < MAX_DIVERGENCE_LIMIT,
        "Velocities should remain nearly divergence-free after long run (max div = {max_divergence})"
    );

    // Verify: all fields remain finite
    assert!(
        grid.cells.iter().all(|cell| cell.rho_o2.is_finite()
            && cell.rho_n2.is_finite()
            && cell.rho_co2.is_finite()),
        "All densities should remain finite"
    );
}

// Helper functions

fn simulate_step(mut atmosphere: ResMut<AtmosphereGrid>, collision_map: Res<TileCollisionMap>) {
    // Run one full compressible step (advection + diffusion + compression)
    fluid_dynamics_research_sim::atmosphere::steps::advection::advection_step(
        &mut atmosphere,
        &collision_map,
        DT,
    );
    fluid_dynamics_research_sim::atmosphere::steps::diffusion::diffusion_step(
        &mut atmosphere,
        &collision_map,
        DT,
    );
    fluid_dynamics_research_sim::atmosphere::steps::compression::compression_heating_step(
        &mut atmosphere,
        &collision_map,
        DT,
    );
    update_grid_pressures(&mut atmosphere);
}

fn compute_co2_center_of_mass(grid: &AtmosphereGrid) -> (f32, f32) {
    let mut total_co2 = 0.0;
    let mut com_x = 0.0;
    let mut com_y = 0.0;

    for y in 0..grid.height {
        for x in 0..grid.width {
            let idx = grid.index(x, y);
            let co2 = grid.cells[idx].rho_co2;
            total_co2 += co2;
            com_x += co2 * x as f32;
            com_y += co2 * y as f32;
        }
    }

    if total_co2 > 0.0 {
        (com_x / total_co2, com_y / total_co2)
    } else {
        (0.0, 0.0)
    }
}

fn compute_max_divergence(grid: &AtmosphereGrid, collision_map: &TileCollisionMap) -> f32 {
    let dx = grid.tile_size_physical;
    let mut max_div: f32 = 0.0;

    for y in 0..grid.height {
        for x in 0..grid.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = grid.index(x, y);
            let cell = &grid.cells[idx];

            let u_right = if x < grid.width - 1 && !collision_map.is_blocked(x + 1, y) {
                grid.cells[grid.index(x + 1, y)].u
            } else {
                cell.u
            };

            let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                grid.cells[grid.index(x - 1, y)].u
            } else {
                cell.u
            };

            let v_up = if y < grid.height - 1 && !collision_map.is_blocked(x, y + 1) {
                grid.cells[grid.index(x, y + 1)].v
            } else {
                cell.v
            };

            let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                grid.cells[grid.index(x, y - 1)].v
            } else {
                cell.v
            };

            let du_dx = (u_right - u_left) / (2.0 * dx);
            let dv_dy = (v_up - v_down) / (2.0 * dx);
            let div = (du_dx + dv_dy).abs();

            max_div = max_div.max(div);
        }
    }

    max_div
}

fn compute_total_mass(grid: &AtmosphereGrid) -> f32 {
    let tile_volume = grid.tile_size_physical * grid.tile_size_physical * constants::ROOM_HEIGHT;
    grid.cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum()
}

fn compute_total_momentum(grid: &AtmosphereGrid) -> (f32, f32) {
    let tile_volume = grid.tile_size_physical * grid.tile_size_physical * constants::ROOM_HEIGHT;
    let mut mom_x = 0.0;
    let mut mom_y = 0.0;

    for cell in grid.cells.iter() {
        let mass = cell.total_density() * tile_volume;
        mom_x += mass * cell.u;
        mom_y += mass * cell.v;
    }

    (mom_x, mom_y)
}
