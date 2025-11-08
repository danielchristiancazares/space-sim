use super::debug::debug_visualization;
use super::grid::AtmosphereGrid;
use super::monitoring::{
    check_mass_conservation, monitor_divergence, DivergenceTracker, MassTracker,
};
use super::simulation::{simulate_atmosphere, update_grid_pressures};
use super::sources::{life_support_generation, life_support_mixing, player_respiration};
use crate::tilemap::{TileCollisionMap, TilemapInitSet};
use bevy::prelude::*;

pub struct AtmospherePlugin;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub enum AtmosphereSimSet {
    Main,
}

impl Plugin for AtmospherePlugin {
    fn build(&self, app: &mut App) {
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
                    check_mass_conservation,
                    monitor_divergence,
                    debug_visualization,
                )
                    .chain()
                    .in_set(AtmosphereSimSet::Main),
            );
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

    // Initialize mass tracker with zero mass (starting in vacuum)
    commands.insert_resource(MassTracker::default());

    // Initialize divergence tracker
    commands.insert_resource(DivergenceTracker::default());

    info!(
        "Atmospheric simulation initialized: {}x{} grid (starting in vacuum)",
        collision_map.width, collision_map.height
    );
}

fn update_pressure_from_state(mut atmosphere: ResMut<AtmosphereGrid>) {
    update_grid_pressures(&mut atmosphere);
}
