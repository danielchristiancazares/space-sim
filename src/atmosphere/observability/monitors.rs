use super::metrics::SimulationDiagnostics;
use super::ObservabilityConfig;
use crate::atmosphere::constants;
use bevy::prelude::*;

#[allow(dead_code)]
#[derive(Message, Debug, Clone)]
pub enum SimulationViolationEvent {
    NaNInf {
        count: u32,
    },
    NegativeDensity {
        count: u32,
    },
    Supersonic {
        count: u32,
        max_velocity: f32,
    },
    TemperatureViolation {
        count: u32,
        min: f32,
        max: f32,
    },
    CFLCrossed {
        value: f32,
    },
    DivergenceGrowth {
        prev: f32,
        current: f32,
        factor: f32,
    },
    DivergenceHigh {
        max: f32,
        avg: f32,
    },
    MassConservation {
        abs: f32,
        relative: f32,
        prev: f32,
        current: f32,
        expected_delta: f32,
    },
    MomentumChange {
        dx: f32,
        dy: f32,
        magnitude: f32,
    },
}

/// Interval-gated emission of violation warnings and events.
pub fn emit_simulation_violations(
    diagnostics: Res<SimulationDiagnostics>,
    config: Res<ObservabilityConfig>,
    time: Res<Time>,
    mut events: Option<ResMut<Messages<SimulationViolationEvent>>>,
) {
    let sim_time = diagnostics.sim_time;
    let interval = config.violation_check_interval;
    if interval <= 0.0 {
        return;
    }
    let time_mod = sim_time % interval;
    if time_mod >= time.delta_secs() {
        return;
    }

    // Log with existing helper
    diagnostics.check_violations();

    if !config.enable_violation_events {
        return;
    }

    let Some(ref mut events) = events else {
        return;
    };

    if diagnostics.nan_inf_count > 0 {
        events.write(SimulationViolationEvent::NaNInf {
            count: diagnostics.nan_inf_count,
        });
    }
    if diagnostics.negative_density_count > 0 {
        events.write(SimulationViolationEvent::NegativeDensity {
            count: diagnostics.negative_density_count,
        });
    }
    if diagnostics.supersonic_count > 0 {
        events.write(SimulationViolationEvent::Supersonic {
            count: diagnostics.supersonic_count,
            max_velocity: diagnostics.max_velocity,
        });
    }
    if diagnostics.temp_violation_count > 0 {
        events.write(SimulationViolationEvent::TemperatureViolation {
            count: diagnostics.temp_violation_count,
            min: diagnostics.min_temperature,
            max: diagnostics.max_temperature,
        });
    }
    if diagnostics.max_cfl > 1.0 {
        events.write(SimulationViolationEvent::CFLCrossed {
            value: diagnostics.max_cfl,
        });
    }
    // Divergence growth warning already logged in metrics; emit event as well
    if diagnostics.div_growth_factor > 2.0 && diagnostics.prev_max_div > 1e-6 {
        events.write(SimulationViolationEvent::DivergenceGrowth {
            prev: diagnostics.prev_max_div,
            current: diagnostics.max_div_u,
            factor: diagnostics.div_growth_factor,
        });
    }
}

/// Periodic detailed logging of diagnostics (1 Hz by default).
pub fn log_detailed_diagnostics(
    diagnostics: Res<SimulationDiagnostics>,
    config: Res<ObservabilityConfig>,
    time: Res<Time>,
) {
    let sim_time = diagnostics.sim_time;
    let interval = config.detailed_log_interval.max(0.0);
    if interval <= 0.0 {
        return;
    }
    let time_mod = sim_time % interval;
    if time_mod < time.delta_secs() {
        diagnostics.log_detailed_stats();
    }
}

/// Mass conservation check using SimulationDiagnostics values; logs and emits event.
pub fn check_mass_conservation(
    mut diagnostics: ResMut<SimulationDiagnostics>,
    mut time_since: Local<f32>,
    time: Res<Time>,
    mut events: Option<ResMut<Messages<SimulationViolationEvent>>>,
) {
    *time_since += time.delta_secs();
    if *time_since < constants::MASS_CHECK_INTERVAL {
        return;
    }

    let total_mass = diagnostics.total_mass();
    if diagnostics.last_total_mass > f32::EPSILON {
        let actual_delta = total_mass - diagnostics.last_total_mass;
        let error = (actual_delta - diagnostics.expected_mass_delta).abs();
        let relative_error = if diagnostics.last_total_mass > f32::EPSILON {
            error / diagnostics.last_total_mass
        } else {
            0.0
        };

        if relative_error > constants::MASS_TOLERANCE {
            warn!(
                "Mass conservation warning: Total mass changed by {:.6} kg (expected {:.6} kg). Relative error: {:.2}%",
                actual_delta,
                diagnostics.expected_mass_delta,
                relative_error * 100.0
            );
            if let Some(ref mut events) = events {
                events.write(SimulationViolationEvent::MassConservation {
                    abs: error,
                    relative: relative_error,
                    prev: diagnostics.last_total_mass,
                    current: total_mass,
                    expected_delta: diagnostics.expected_mass_delta,
                });
            }
        } else {
            debug!(
                "Mass conservation check: Total mass = {:.3} kg (delta: {:.6} kg, expected: {:.6} kg)",
                total_mass, actual_delta, diagnostics.expected_mass_delta
            );
        }
    }

    diagnostics.last_total_mass = total_mass;
    diagnostics.expected_mass_delta = 0.0;
    *time_since = 0.0;
}

/// Divergence monitor using precomputed diagnostics.
pub fn monitor_divergence(
    diagnostics: Res<SimulationDiagnostics>,
    mut time_since: Local<f32>,
    time: Res<Time>,
    mut events: Option<ResMut<Messages<SimulationViolationEvent>>>,
) {
    *time_since += time.delta_secs();
    if *time_since < constants::DIVERGENCE_CHECK_INTERVAL {
        return;
    }

    debug!(
        "Divergence check: max = {:.6} s⁻¹, avg = {:.6} s⁻¹",
        diagnostics.max_divergence, diagnostics.avg_divergence
    );

    if diagnostics.max_divergence > 1.0 {
        warn!(
            "High velocity divergence: {:.3} s⁻¹ (compressible transient or CFL issue if sustained)",
            diagnostics.max_divergence
        );
        if let Some(ref mut events) = events {
            events.write(SimulationViolationEvent::DivergenceHigh {
                max: diagnostics.max_divergence,
                avg: diagnostics.avg_divergence,
            });
        }
    }

    *time_since = 0.0;
}

/// Momentum drift monitor using precomputed diagnostics.
pub fn check_momentum_conservation(
    mut diagnostics: ResMut<SimulationDiagnostics>,
    mut time_since: Local<f32>,
    time: Res<Time>,
    mut events: Option<ResMut<Messages<SimulationViolationEvent>>>,
) {
    *time_since += time.delta_secs();
    if *time_since < constants::MASS_CHECK_INTERVAL {
        return;
    }

    let total_momentum_x = diagnostics.momentum_x;
    let total_momentum_y = diagnostics.momentum_y;

    if diagnostics.last_momentum_x.abs() > f32::EPSILON
        || diagnostics.last_momentum_y.abs() > f32::EPSILON
    {
        let delta_x = total_momentum_x - diagnostics.last_momentum_x;
        let delta_y = total_momentum_y - diagnostics.last_momentum_y;
        let delta_magnitude = (delta_x * delta_x + delta_y * delta_y).sqrt();

        if delta_magnitude > 10.0 {
            debug!(
                "Momentum change: Δp = ({:.3}, {:.3}) kg·m/s (magnitude: {:.3})",
                delta_x, delta_y, delta_magnitude
            );
            if let Some(ref mut events) = events {
                events.write(SimulationViolationEvent::MomentumChange {
                    dx: delta_x,
                    dy: delta_y,
                    magnitude: delta_magnitude,
                });
            }
        }
    }

    diagnostics.last_momentum_x = total_momentum_x;
    diagnostics.last_momentum_y = total_momentum_y;
    *time_since = 0.0;
}
