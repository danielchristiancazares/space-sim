pub mod constants;

mod grid;
mod plugin;
mod simulation;

pub mod steps {
    pub mod advection;
    pub mod diffusion;
    pub mod pressure_flux;
    pub mod pressure_projection;
}

mod debug;
mod monitoring;
mod sources;

#[allow(unused_imports)]
pub use grid::{AtmosphereCell, AtmosphereGrid};
#[allow(unused_imports)]
pub use monitoring::{DivergenceTracker, MassTracker};
pub use plugin::{AtmospherePlugin, AtmosphereSimSet};
#[allow(unused_imports)]
pub use simulation::update_grid_pressures;
