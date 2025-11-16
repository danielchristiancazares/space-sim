use super::grid::AtmosphereGrid;
use super::observability::{AtmosphereObservabilityPlugin, CsvExporter, SimulationDiagnostics};
// Monitoring systems are provided by the observability plugin
use super::simulation::{simulate_atmosphere, update_grid_pressures};
use super::sources::{
    life_support_generation, life_support_mixing, player_respiration, BreathabilityWarningTracker,
};
use crate::tilemap::{TileCollisionMap, TilemapInitSet};
use bevy::prelude::*;

pub struct AtmospherePlugin;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub enum AtmosphereSimSet {
    Main,
}

impl Plugin for AtmospherePlugin {
    fn build(&self, app: &mut App) {
        // Research mode: Run simulation in Update schedule at full CPU speed
        // (no fixed timestep rate limiting - simulation runs as fast as possible)
        app.configure_sets(Update, AtmosphereSimSet::Main);
        app.add_systems(Startup, initialize_atmosphere.after(TilemapInitSet))
            .add_systems(
                Update,
                (
                    life_support_generation,
                    life_support_mixing,
                    player_respiration,
                    update_pressure_from_state,
                    simulate_atmosphere,
                )
                    .chain()
                    .in_set(AtmosphereSimSet::Main),
            );
        // Add observability wiring (monitors, visualization, pressure logger)
        app.add_plugins(AtmosphereObservabilityPlugin);
    }
}

fn initialize_atmosphere(mut commands: Commands, collision_map: Res<TileCollisionMap>) {
    let mut grid = AtmosphereGrid::new(
        collision_map.width,
        collision_map.height,
        collision_map.tile_size.x,
        1.0, // Treat each tile as a 1m × 1m column of air
    );

    // Start with vacuum - life support will fill the room
    grid.initialize_vacuum();

    commands.insert_resource(grid);

    // Initialize breathability warning tracker
    commands.insert_resource(BreathabilityWarningTracker::default());

    // Initialize diagnostics and CSV exporter (writes to .debug/diagnostics/ every 5 seconds)
    commands.insert_resource(SimulationDiagnostics::default());
    commands.insert_resource(CsvExporter::default());

    info!(
        "Atmospheric simulation initialized: {}x{} grid (starting in vacuum)",
        collision_map.width, collision_map.height
    );
    info!("Diagnostics CSV export: .debug/diagnostics/ (5 second intervals)");
}

fn update_pressure_from_state(mut atmosphere: ResMut<AtmosphereGrid>) {
    update_grid_pressures(&mut atmosphere);
}
