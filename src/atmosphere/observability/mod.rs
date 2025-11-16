pub mod export;
pub mod metrics;
pub mod monitors;
pub mod pressure_logger;
pub mod viz;

use bevy::prelude::*;

// Re-exports for existing call sites
pub use export::{estimate_file_size, CsvExporter};
pub use metrics::SimulationDiagnostics;
pub use pressure_logger::PressureLoggerPlugin;
pub use viz::debug_visualization;

/// Toggles and intervals for observability features.
#[derive(Resource, Clone)]
pub struct ObservabilityConfig {
    pub enable_visualization: bool,
    pub enable_pressure_logger: bool,
    pub enable_csv_export: bool,
    pub enable_violation_events: bool,
    pub detailed_log_interval: f32,
    pub violation_check_interval: f32,
}

impl Default for ObservabilityConfig {
    fn default() -> Self {
        Self {
            enable_visualization: true,
            enable_pressure_logger: true,
            enable_csv_export: true,
            enable_violation_events: true,
            detailed_log_interval: 1.0,
            violation_check_interval: 5.0,
        }
    }
}

/// Plugin that wires diagnostics/monitoring/visualization systems.
pub struct AtmosphereObservabilityPlugin;

impl Plugin for AtmosphereObservabilityPlugin {
    fn build(&self, app: &mut App) {
        use crate::atmosphere::plugin::AtmosphereSimSet;

        // Default config unless provided by caller
        if app.world().get_resource::<ObservabilityConfig>().is_none() {
            app.insert_resource(ObservabilityConfig::default());
        }

        app.add_message::<monitors::SimulationViolationEvent>();

        // Snapshot config to decide which systems to register
        let config = app
            .world()
            .get_resource::<ObservabilityConfig>()
            .cloned()
            .unwrap_or_default();

        // Always wire core monitors
        app.add_systems(
            Update,
            (
                monitors::emit_simulation_violations,
                monitors::log_detailed_diagnostics,
                monitors::check_mass_conservation,
                monitors::check_momentum_conservation,
                monitors::monitor_divergence,
            )
                .chain()
                .in_set(AtmosphereSimSet::Main),
        );

        // Visualization (toggle at registration time)
        if config.enable_visualization {
            app.add_systems(
                Update,
                viz::debug_visualization.in_set(AtmosphereSimSet::Main),
            );
        }

        // Eagerly read other toggles to avoid dead_code warnings
        let _ = config.enable_csv_export;

        // Pressure logger is kept as a separate Bevy plugin
        if config.enable_pressure_logger {
            app.add_plugins(PressureLoggerPlugin);
        }
    }
}
