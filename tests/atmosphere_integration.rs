use bevy::log::{Level, LogPlugin};
use bevy::prelude::*;
use bevy::time::TimeUpdateStrategy;
use fluid_dynamics_research_sim::atmosphere::steps::diffusion::diffusion_step;
use fluid_dynamics_research_sim::{
    atmosphere::{
        constants, life_support_generation, life_support_mixing, simulate_atmosphere,
        update_grid_pressures, AtmosphereGrid, CsvExporter, SimulationDiagnostics,
    },
    tilemap::{LifeSupportTiles, TileCollisionMap},
};
use std::time::Duration;

const GRID_WIDTH: u32 = 7;
const GRID_HEIGHT: u32 = 7;
const DT: f32 = 0.05;
const STEPS: usize = 120;
const MASS_TOLERANCE: f32 = 0.05;

fn assert_all_finite(label: &str, grid: &AtmosphereGrid) {
    for (i, cell) in grid.cells.iter().enumerate() {
        if !cell.u.is_finite()
            || !cell.v.is_finite()
            || !cell.temperature.is_finite()
            || !cell.rho_o2.is_finite()
            || !cell.rho_n2.is_finite()
            || !cell.rho_co2.is_finite()
        {
            let x = i as u32 % GRID_WIDTH;
            let y = i as u32 / GRID_WIDTH;
            panic!(
                "{}: NaN/Inf at cell ({}, {}): rho_o2={}, rho_n2={}, rho_co2={}, u={}, v={}, T={}",
                label,
                x,
                y,
                cell.rho_o2,
                cell.rho_n2,
                cell.rho_co2,
                cell.u,
                cell.v,
                cell.temperature
            );
        }
    }
}

#[test]
fn life_support_fill_remains_stable() {
    let mut app = App::new();
    app.add_plugins(MinimalPlugins);
    app.add_plugins(LogPlugin {
        level: Level::DEBUG,
        filter: "fluid_dynamics_research_sim=debug".into(),
        ..default()
    });

    let collision_map = TileCollisionMap::empty(GRID_WIDTH, GRID_HEIGHT, Vec2::ONE, Vec2::ZERO);
    app.insert_resource(collision_map);
    app.insert_resource(LifeSupportTiles {
        positions: vec![UVec2::new(GRID_WIDTH / 2, GRID_HEIGHT / 2)],
    });

    let mut grid = AtmosphereGrid::new(GRID_WIDTH, GRID_HEIGHT, 1.0, 1.0);
    grid.initialize_vacuum();
    app.insert_resource(grid);

    app.insert_resource(SimulationDiagnostics::default());
    app.insert_resource(CsvExporter::default());
    app.world_mut()
        .insert_resource(TimeUpdateStrategy::ManualDuration(Duration::from_secs_f32(
            DT,
        )));

    app.add_systems(
        Update,
        (
            life_support_generation,
            life_support_mixing,
            sync_pressures,
            simulate_atmosphere,
            sync_pressures_post,
        )
            .chain(),
    );

    for step in 0..STEPS {
        app.update();

        {
            let grid = app.world().resource::<AtmosphereGrid>();
            assert_all_finite(&format!("after step {}", step + 1), &grid);
        }
    }

    let grid = app.world().resource::<AtmosphereGrid>();

    assert!(grid.cells.iter().all(|cell| cell.rho_o2.is_finite()
        && cell.rho_n2.is_finite()
        && cell.rho_co2.is_finite()));
    assert!(grid
        .cells
        .iter()
        .all(|cell| cell.u.is_finite() && cell.v.is_finite() && cell.temperature.is_finite()));
    assert!(grid.cells.iter().all(|cell| cell.pressure.is_finite()));

    let tile_volume = grid.tile_size_physical * grid.tile_size_physical * constants::ROOM_HEIGHT;
    let total_mass: f32 = grid
        .cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum();
    let expected_mass =
        (constants::LIFE_SUPPORT_O2_RATE + constants::LIFE_SUPPORT_N2_RATE) * DT * STEPS as f32;
    let tolerance = expected_mass * MASS_TOLERANCE;
    debug!(
        "Integration final mass {:.3} kg (expected {:.3} ± {:.3})",
        total_mass, expected_mass, tolerance
    );
    assert!(
        (total_mass - expected_mass).abs() <= tolerance,
        "Total mass {total_mass} deviated from expected {expected_mass}"
    );

    let max_pressure = grid
        .cells
        .iter()
        .map(|cell| cell.pressure)
        .fold(0.0_f32, f32::max);
    assert!(max_pressure > 0.0, "Atmosphere failed to build pressure");

    let center_idx = grid.index(GRID_WIDTH / 2, GRID_HEIGHT / 2);
    let center_pressure = grid.cells[center_idx].pressure;
    debug!(
        "Integration pressure summary: center {:.1} Pa, max {:.1} Pa",
        center_pressure, max_pressure
    );
    assert!(
        center_pressure >= max_pressure - 1.0,
        "Center pressure {center_pressure} not near max {max_pressure}"
    );
}

fn sync_pressures(mut atmosphere: ResMut<AtmosphereGrid>) {
    update_grid_pressures(&mut atmosphere);
}

fn sync_pressures_post(mut atmosphere: ResMut<AtmosphereGrid>) {
    update_grid_pressures(&mut atmosphere);
}

#[test]
fn diffusion_into_vacuum_equalizes_mass() {
    let mut grid = AtmosphereGrid::new(2, 1, 1.0, 1.0);
    let collision_map = TileCollisionMap::empty(2, 1, Vec2::ONE, Vec2::ZERO);

    grid.initialize_vacuum();

    let left = grid.index(0, 0);
    let right = grid.index(1, 0);
    let tile_volume = grid.tile_size_physical * grid.tile_size_physical * constants::ROOM_HEIGHT;

    // Seed the left cell with Earth-like oxygen density, right cell stays vacuum.
    grid.cells[left].rho_o2 = constants::EARTH_RHO_O2;
    grid.cells[left].update_pressure();
    grid.cells[right].update_pressure();

    let initial_mass = grid.cells[left].total_density() * tile_volume;
    assert!(initial_mass > 0.0);

    // Use a single diffusion step with a stable timestep (alpha = 0.05 < 0.2).
    let dt = 100.0;
    diffusion_step(&mut grid, &collision_map, dt);
    update_grid_pressures(&mut grid);

    let alpha = (constants::GAS_DIFFUSIVITY * dt
        / (grid.tile_size_physical * grid.tile_size_physical))
        .min(0.2);
    // With zero-gradient (Neumann) boundaries, total mass is conserved.
    // For a 2-cell domain, the discrete update gives:
    // left_mass' = left_mass * (1 - alpha)
    // right_mass' = left_mass * alpha
    let expected_left = initial_mass * (1.0 - alpha);
    let expected_right = initial_mass * alpha;
    let expected_total = initial_mass;
    let tolerance = initial_mass * 0.02; // Allow 2% numerical wiggle room

    let total_mass_after: f32 = grid
        .cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum();
    assert!(
        (total_mass_after - expected_total).abs() <= tolerance,
        "Total mass {total_mass_after} deviated from expected {expected_total}"
    );

    let left_mass = grid.cells[left].total_density() * tile_volume;
    let right_mass = grid.cells[right].total_density() * tile_volume;
    assert!(
        (left_mass - expected_left).abs() <= tolerance,
        "Left mass {left_mass} did not match expected {expected_left}"
    );
    assert!(
        (right_mass - expected_right).abs() <= tolerance,
        "Right mass {right_mass} did not match expected {expected_right}"
    );

    // Ensure the vacuum cell actually gained mass.
    assert!(
        right_mass > 0.0,
        "Right cell should gain mass from diffusion into vacuum"
    );
}

// Legacy flux/projection comparison test removed: pressure-driven flux path no longer used.
