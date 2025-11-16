use crate::atmosphere::constants;
use crate::atmosphere::grid::AtmosphereGrid;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::{info, warn, Resource};

/// Lightweight diagnostics tracking fundamental simulation correctness.
/// Focused on "sanity checks" rather than high-precision error tracking.
#[derive(Resource)]
pub struct SimulationDiagnostics {
    /// Current simulation frame number
    pub frame: u64,
    /// Total accumulated simulation time [s]
    pub sim_time: f32,

    // Conservation quantities
    /// Total O2 mass [kg]
    pub mass_o2: f32,
    /// Total N2 mass [kg]
    pub mass_n2: f32,
    /// Total CO2 mass [kg]
    pub mass_co2: f32,
    /// Total momentum in x-direction [kg·m/s]
    pub momentum_x: f32,
    /// Total momentum in y-direction [kg·m/s]
    pub momentum_y: f32,
    /// Total kinetic energy [J]
    pub kinetic_energy: f32,

    // Physical extremes (for sanity checks)
    /// Maximum velocity magnitude [m/s]
    pub max_velocity: f32,
    /// Maximum temperature [K]
    pub max_temperature: f32,
    /// Minimum temperature [K]
    pub min_temperature: f32,
    /// Minimum total density [kg/m³]
    pub min_density: f32,
    /// Maximum total density [kg/m³]
    pub max_density: f32,

    // Solver quality metrics
    /// Maximum velocity divergence [s⁻¹]
    pub max_divergence: f32,
    /// Average velocity divergence [s⁻¹]
    pub avg_divergence: f32,
    /// Snapshot divergence recorded before any optional projection instrumentation [s⁻¹]
    pub divergence_before_projection: f32,
    /// Iteration count reported by optional/legacy pressure projection instrumentation
    pub pressure_iterations: usize,

    // Timestep tracking
    /// Current timestep [s]
    pub timestep: f32,
    /// Maximum CFL number encountered
    pub max_cfl: f32,

    // Violation counts (for warnings)
    pub nan_inf_count: u32,
    pub negative_density_count: u32,
    pub supersonic_count: u32,
    pub temp_violation_count: u32,

    // Conservation tracking (research-grade accuracy)
    /// Mass from previous step [kg] (to compute Δmass)
    pub prev_total_mass: f32,
    /// Momentum magnitude from previous step [kg·m/s] (to compute Δmomentum)
    pub prev_momentum_mag: f32,
    /// Current step mass error [kg]
    pub mass_error: f32,
    /// Current step momentum error [kg·m/s]
    pub momentum_error: f32,
    /// Maximum mass error observed [kg]
    pub max_mass_error: f32,
    /// Maximum momentum error observed [kg·m/s]
    pub max_momentum_error: f32,
    /// Total accumulated mass error [kg]
    pub cumulative_mass_error: f32,
    /// Total accumulated momentum error [kg·m/s]
    pub cumulative_momentum_error: f32,

    // Detailed statistics (temperature spike investigation)
    /// Average temperature across domain [K]
    pub avg_temperature: f32,
    /// Maximum velocity divergence [s⁻¹]
    pub max_div_u: f32,
    /// Minimum velocity divergence (most negative = strongest compression) [s⁻¹]
    pub min_div_u: f32,
    /// Average velocity divergence [s⁻¹]
    pub avg_div_u: f32,
    /// Average compression heating rate [K/s]
    pub compression_heating_rate: f32,

    // Divergence growth tracking (instability detection)
    /// Previous maximum divergence [s⁻¹] (for growth rate calculation)
    pub prev_max_div: f32,
    /// Divergence growth factor: max_div_current / max_div_prev
    pub div_growth_factor: f32,
}

impl Default for SimulationDiagnostics {
    fn default() -> Self {
        Self {
            frame: 0,
            sim_time: 0.0,
            mass_o2: 0.0,
            mass_n2: 0.0,
            mass_co2: 0.0,
            momentum_x: 0.0,
            momentum_y: 0.0,
            kinetic_energy: 0.0,
            max_velocity: 0.0,
            max_temperature: 0.0,
            min_temperature: f32::INFINITY,
            min_density: f32::INFINITY,
            max_density: 0.0,
            max_divergence: 0.0,
            avg_divergence: 0.0,
            divergence_before_projection: 0.0,
            pressure_iterations: 0,
            timestep: 0.0,
            max_cfl: 0.0,
            nan_inf_count: 0,
            negative_density_count: 0,
            supersonic_count: 0,
            temp_violation_count: 0,
            prev_total_mass: 0.0,
            prev_momentum_mag: 0.0,
            mass_error: 0.0,
            momentum_error: 0.0,
            max_mass_error: 0.0,
            max_momentum_error: 0.0,
            cumulative_mass_error: 0.0,
            cumulative_momentum_error: 0.0,
            avg_temperature: 0.0,
            max_div_u: 0.0,
            min_div_u: 0.0,
            avg_div_u: 0.0,
            compression_heating_rate: 0.0,
            prev_max_div: 0.0,
            div_growth_factor: 1.0,
        }
    }
}

impl SimulationDiagnostics {
    /// Compute all metrics from current atmosphere state.
    pub fn update(
        &mut self,
        atmosphere: &AtmosphereGrid,
        collision_map: &TileCollisionMap,
        dt: f32,
    ) {
        self.frame += 1;
        self.sim_time += dt;
        self.timestep = dt;

        // Reset aggregates
        self.mass_o2 = 0.0;
        self.mass_n2 = 0.0;
        self.mass_co2 = 0.0;
        self.momentum_x = 0.0;
        self.momentum_y = 0.0;
        self.kinetic_energy = 0.0;
        self.max_velocity = 0.0;
        self.max_temperature = 0.0;
        self.min_temperature = f32::INFINITY;
        self.min_density = f32::INFINITY;
        self.max_density = 0.0;
        self.nan_inf_count = 0;
        self.negative_density_count = 0;
        self.supersonic_count = 0;
        self.temp_violation_count = 0;

        // Reset detailed stats
        let mut sum_temperature = 0.0;
        self.max_div_u = f32::NEG_INFINITY;
        self.min_div_u = f32::INFINITY;
        let mut stat_count = 0;

        let tile_volume =
            atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

        // Speed of sound in air at room temp [m/s]
        const SPEED_OF_SOUND: f32 = 343.0;

        for y in 0..atmosphere.height {
            for x in 0..atmosphere.width {
                if collision_map.is_blocked(x, y) {
                    continue;
                }

                let idx = atmosphere.index(x, y);
                let cell = &atmosphere.cells[idx];

                // Check for NaN/Inf
                if !cell.rho_o2.is_finite()
                    || !cell.rho_n2.is_finite()
                    || !cell.rho_co2.is_finite()
                    || !cell.u.is_finite()
                    || !cell.v.is_finite()
                    || !cell.temperature.is_finite()
                {
                    self.nan_inf_count += 1;
                    continue; // Skip this cell's contributions
                }

                let rho = cell.total_density();
                let mass = rho * tile_volume;

                // Check for negative densities
                if cell.rho_o2 < 0.0 || cell.rho_n2 < 0.0 || cell.rho_co2 < 0.0 {
                    self.negative_density_count += 1;
                }

                // Conservation quantities
                self.mass_o2 += cell.rho_o2 * tile_volume;
                self.mass_n2 += cell.rho_n2 * tile_volume;
                self.mass_co2 += cell.rho_co2 * tile_volume;
                self.momentum_x += mass * cell.u;
                self.momentum_y += mass * cell.v;

                let velocity_sq = cell.u * cell.u + cell.v * cell.v;
                let velocity = velocity_sq.sqrt();
                self.kinetic_energy += 0.5 * mass * velocity_sq;

                // Extremes
                self.max_velocity = self.max_velocity.max(velocity);
                self.max_temperature = self.max_temperature.max(cell.temperature);
                self.min_temperature = self.min_temperature.min(cell.temperature);
                self.min_density = self.min_density.min(rho);
                self.max_density = self.max_density.max(rho);

                // Accumulate for averages
                sum_temperature += cell.temperature;
                stat_count += 1;

                // Sanity checks
                if velocity > SPEED_OF_SOUND {
                    self.supersonic_count += 1;
                }

                if cell.temperature < constants::T_CMB || cell.temperature > constants::T_PHYS_MAX {
                    self.temp_violation_count += 1;
                }
            }
        }

        // Compute CFL number
        if self.max_velocity > 0.0 {
            self.max_cfl = (self.max_velocity * dt) / atmosphere.tile_size_physical;
        } else {
            self.max_cfl = 0.0;
        }

        // Compute conservation errors (research-grade accuracy tracking)
        let current_mass = self.total_mass();
        let current_momentum = self.momentum_magnitude();

        if self.frame > 1 {
            // Compute change in mass and momentum (should be near zero for conservation)
            self.mass_error = (current_mass - self.prev_total_mass).abs();
            self.momentum_error = (current_momentum - self.prev_momentum_mag).abs();

            // Track maximum errors observed
            self.max_mass_error = self.max_mass_error.max(self.mass_error);
            self.max_momentum_error = self.max_momentum_error.max(self.momentum_error);

            // Accumulate total error over simulation
            self.cumulative_mass_error += self.mass_error;
            self.cumulative_momentum_error += self.momentum_error;

            // Log significant conservation violations
            let relative_mass_error = if self.prev_total_mass > 0.0 {
                self.mass_error / self.prev_total_mass
            } else {
                0.0
            };

            if relative_mass_error > 1e-6 {
                bevy::log::debug!(
                    "Frame {}: Mass conservation error: Δ={:.6e} kg ({:.3e} relative), \
                     prev={:.6} kg, current={:.6} kg",
                    self.frame,
                    self.mass_error,
                    relative_mass_error,
                    self.prev_total_mass,
                    current_mass
                );
            }
        }

        // Store current values for next frame comparison
        self.prev_total_mass = current_mass;
        self.prev_momentum_mag = current_momentum;

        // Finalize averages
        if stat_count > 0 {
            self.avg_temperature = sum_temperature / stat_count as f32;
        }

        // Compute divergence statistics
        self.compute_divergence_statistics(atmosphere, collision_map);

        // Estimate compression heating rate: -(γ-1) × T_avg × div_u_avg
        self.compression_heating_rate = -(constants::GAMMA - 1.0) * self.avg_temperature * self.avg_div_u;

        // Old divergence computation (keep for compatibility)
        self.compute_divergence(atmosphere, collision_map);
    }

    /// Compute velocity divergence field statistics.
    fn compute_divergence(
        &mut self,
        atmosphere: &AtmosphereGrid,
        collision_map: &TileCollisionMap,
    ) {
        let dx = atmosphere.tile_size_physical;
        let mut max_div: f32 = 0.0;
        let mut total_div = 0.0;
        let mut count = 0;

        for y in 0..atmosphere.height {
            for x in 0..atmosphere.width {
                if collision_map.is_blocked(x, y) {
                    continue;
                }

                let idx = atmosphere.index(x, y);
                let cell = &atmosphere.cells[idx];

                // Central difference for divergence
                let u_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                    atmosphere.cells[atmosphere.index(x + 1, y)].u
                } else {
                    cell.u
                };

                let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                    atmosphere.cells[atmosphere.index(x - 1, y)].u
                } else {
                    cell.u
                };

                let v_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                    atmosphere.cells[atmosphere.index(x, y + 1)].v
                } else {
                    cell.v
                };

                let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                    atmosphere.cells[atmosphere.index(x, y - 1)].v
                } else {
                    cell.v
                };

                let du_dx = (u_right - u_left) / (2.0 * dx);
                let dv_dy = (v_up - v_down) / (2.0 * dx);
                let div = du_dx + dv_dy;

                max_div = max_div.max(div.abs());
                total_div += div.abs();
                count += 1;
            }
        }

        self.max_divergence = max_div;
        self.avg_divergence = if count > 0 {
            total_div / count as f32
        } else {
            0.0
        };
    }

    /// Compute detailed divergence statistics (min, max, avg with sign).
    fn compute_divergence_statistics(
        &mut self,
        atmosphere: &AtmosphereGrid,
        collision_map: &TileCollisionMap,
    ) {
        let dx = atmosphere.tile_size_physical;
        let mut sum_div = 0.0;
        let mut count = 0;
        let mut max_div = f32::NEG_INFINITY;
        let mut min_div = f32::INFINITY;

        for y in 0..atmosphere.height {
            for x in 0..atmosphere.width {
                if collision_map.is_blocked(x, y) {
                    continue;
                }

                // Compute divergence (with sign)
                let u_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                    atmosphere.cells[atmosphere.index(x + 1, y)].u
                } else {
                    0.0
                };

                let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                    atmosphere.cells[atmosphere.index(x - 1, y)].u
                } else {
                    0.0
                };

                let v_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                    atmosphere.cells[atmosphere.index(x, y + 1)].v
                } else {
                    0.0
                };

                let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                    atmosphere.cells[atmosphere.index(x, y - 1)].v
                } else {
                    0.0
                };

                let du_dx = (u_right - u_left) / (2.0 * dx);
                let dv_dy = (v_up - v_down) / (2.0 * dx);
                let div = du_dx + dv_dy;

                max_div = max_div.max(div);
                min_div = min_div.min(div);
                sum_div += div;
                count += 1;
            }
        }

        // Capture previous divergence BEFORE overwriting for growth tracking
        let current_max_div = max_div;
        let prev_max_div = self.prev_max_div;

        self.max_div_u = max_div;
        self.min_div_u = min_div;
        self.avg_div_u = if count > 0 {
            sum_div / count as f32
        } else {
            0.0
        };

        // ============================================================
        // Divergence growth tracking (instability early warning)
        // ============================================================
        // Tracks how fast divergence is growing frame-to-frame.
        // Rapid growth (> 2×) indicates potential numerical instability.
        self.div_growth_factor = if prev_max_div > 1e-6 {
            current_max_div / prev_max_div
        } else {
            1.0
        };

        self.prev_max_div = current_max_div;

        // Warn if divergence is growing rapidly AND is significant in magnitude
        if self.div_growth_factor > 2.0 && current_max_div > 1e-3 {
            bevy::log::warn!(
                "Divergence growing rapidly: {:.2} → {:.2} s⁻¹ (×{:.2} growth)",
                prev_max_div,
                current_max_div,
                self.div_growth_factor
            );
        }
    }

    /// Record divergence before a (currently optional) pressure projection step.
    #[allow(dead_code)]
    pub fn record_pre_projection_divergence(&mut self, divergence: f32) {
        self.divergence_before_projection = divergence;
    }

    /// Record pressure projection iteration count (legacy instrumentation hook).
    #[allow(dead_code)]
    pub fn record_pressure_iterations(&mut self, iters: usize) {
        self.pressure_iterations = iters;
    }

    /// Check for violations and emit warnings (call this at logging intervals, not every frame).
    pub fn check_violations(&self) {
        if self.nan_inf_count > 0 {
            warn!(
                "Frame {}: NaN/Inf detected in {} cells",
                self.frame, self.nan_inf_count
            );
        }

        if self.negative_density_count > 0 {
            warn!(
                "Frame {}: Negative density in {} cells",
                self.frame, self.negative_density_count
            );
        }

        if self.supersonic_count > 0 {
            warn!(
                "Frame {}: Supersonic velocity (>343 m/s) in {} cells, max = {:.1} m/s",
                self.frame, self.supersonic_count, self.max_velocity
            );
        }

        if self.temp_violation_count > 0 {
            warn!(
                "Frame {}: Temperature out of physical range [2.7 K, 1000.0 K] in {} cells (min={:.1} K, max={:.1} K)",
                self.frame, self.temp_violation_count, self.min_temperature, self.max_temperature
            );
        }

        // Check if pressure projection instrumentation is active and note poor performance.
        if self.divergence_before_projection > 0.0 {
            let reduction = (self.divergence_before_projection - self.max_divergence)
                / self.divergence_before_projection;

            if reduction < 0.5 {
                warn!(
                    "Frame {}: Pressure projection instrumentation reported only {:.1}% divergence reduction. Before: {:.6}, After: {:.6}, Iterations: {}",
                    self.frame,
                    reduction * 100.0,
                    self.divergence_before_projection,
                    self.max_divergence,
                    self.pressure_iterations
                );
            }
        }

        if self.max_cfl > 1.0 {
            warn!(
                "Frame {}: CFL condition violated (CFL = {:.3} > 1.0)",
                self.frame, self.max_cfl
            );
        }
    }

    /// Get total mass for all species combined [kg].
    #[allow(dead_code)]
    pub fn total_mass(&self) -> f32 {
        self.mass_o2 + self.mass_n2 + self.mass_co2
    }

    /// Get total momentum magnitude [kg·m/s].
    #[allow(dead_code)]
    pub fn momentum_magnitude(&self) -> f32 {
        (self.momentum_x * self.momentum_x + self.momentum_y * self.momentum_y).sqrt()
    }

    /// Log detailed simulation statistics (call periodically, not every frame).
    pub fn log_detailed_stats(&self) {
        info!("=== Simulation Stats (t={:.2}s, frame={}) ===", self.sim_time, self.frame);
        info!(
            "  Temperature: avg={:.1} K, min={:.1} K, max={:.1} K",
            self.avg_temperature, self.min_temperature, self.max_temperature
        );
        info!(
            "  Divergence: avg={:.3} s⁻¹, min={:.3} s⁻¹ (compression), max={:.3} s⁻¹ (expansion)",
            self.avg_div_u, self.min_div_u, self.max_div_u
        );
        info!(
            "  Compression heating: avg_rate={:.1} K/s (+ = heating, - = cooling)",
            self.compression_heating_rate
        );
        info!(
            "  Mass: total={:.3} kg (O2={:.3}, N2={:.3}, CO2={:.6})",
            self.total_mass(), self.mass_o2, self.mass_n2, self.mass_co2
        );
        info!(
            "  Conservation: mass_err={:.3e} kg, momentum_err={:.3e} kg·m/s",
            self.mass_error, self.momentum_error
        );
    }
}
