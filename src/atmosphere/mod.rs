#![allow(unused_imports)]

pub mod constants;

mod grid;
pub mod observability;
mod plugin;
mod simulation;

pub mod steps {
    pub mod advection;
    pub mod compression;
    pub mod diffusion;
}

mod sources;

pub use grid::{AtmosphereCell, AtmosphereGrid};
pub use observability::{
    estimate_file_size, CsvExporter, PressureLoggerPlugin, SimulationDiagnostics,
};
pub use plugin::{AtmospherePlugin, AtmosphereSimSet};
pub use simulation::{simulate_atmosphere, update_grid_pressures};
pub use sources::{
    life_support_generation, life_support_mixing, player_respiration, BreathabilityWarningTracker,
};
