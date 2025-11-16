use crate::atmosphere::diagnostics::SimulationDiagnostics;
use bevy::prelude::{info, warn, Resource};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::PathBuf;

/// Manages CSV export of diagnostics data for post-simulation analysis.
/// Writes data at configurable intervals to keep file sizes reasonable.
#[derive(Resource)]
pub struct CsvExporter {
    /// Output directory for CSV files
    output_dir: PathBuf,
    /// File handle for main metrics CSV
    metrics_file: Option<File>,
    /// File handle for conservation errors CSV
    errors_file: Option<File>,
    /// Interval between writes [seconds]
    write_interval: f32,
    /// Time since last write [seconds]
    time_since_write: f32,
    /// Has the header been written?
    headers_written: bool,
    /// Track initial values for error computation
    initial_mass_o2: Option<f32>,
    initial_mass_n2: Option<f32>,
    initial_mass_co2: Option<f32>,
    initial_momentum_x: Option<f32>,
    initial_momentum_y: Option<f32>,
    initial_energy: Option<f32>,
}

impl Default for CsvExporter {
    fn default() -> Self {
        Self::new(".debug/diagnostics", 5.0) // Write every 5 seconds by default
    }
}

impl CsvExporter {
    /// Create a new CSV exporter with specified output directory and write interval.
    ///
    /// # Arguments
    /// * `output_dir` - Directory path for CSV files (created if doesn't exist)
    /// * `write_interval` - Seconds between writes (default 5.0 for ~720 entries/hour)
    pub fn new(output_dir: impl Into<PathBuf>, write_interval: f32) -> Self {
        Self {
            output_dir: output_dir.into(),
            metrics_file: None,
            errors_file: None,
            write_interval,
            time_since_write: 0.0,
            headers_written: false,
            initial_mass_o2: None,
            initial_mass_n2: None,
            initial_mass_co2: None,
            initial_momentum_x: None,
            initial_momentum_y: None,
            initial_energy: None,
        }
    }

    /// Initialize CSV files (creates directory and files with headers).
    pub fn initialize(&mut self) -> io::Result<()> {
        // Create output directory
        fs::create_dir_all(&self.output_dir)?;

        let metrics_path = self.output_dir.join("sim_metrics.csv");
        let errors_path = self.output_dir.join("conservation_errors.csv");

        // Open/create files
        self.metrics_file = Some(
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(&metrics_path)?,
        );

        self.errors_file = Some(
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(&errors_path)?,
        );

        info!(
            "CSV export initialized: metrics={}, errors={}",
            metrics_path.display(),
            errors_path.display()
        );

        Ok(())
    }

    /// Update and potentially write diagnostics data to CSV files.
    pub fn update(&mut self, diagnostics: &SimulationDiagnostics, dt: f32) {
        self.time_since_write += dt;

        // Only write at specified intervals
        if self.time_since_write < self.write_interval {
            return;
        }

        // Initialize files on first write
        if !self.headers_written {
            if let Err(e) = self.initialize() {
                warn!("Failed to initialize CSV export: {}", e);
                return;
            }

            // Write headers
            if let Err(e) = self.write_headers() {
                warn!("Failed to write CSV headers: {}", e);
                return;
            }

            // Capture initial values
            self.initial_mass_o2 = Some(diagnostics.mass_o2);
            self.initial_mass_n2 = Some(diagnostics.mass_n2);
            self.initial_mass_co2 = Some(diagnostics.mass_co2);
            self.initial_momentum_x = Some(diagnostics.momentum_x);
            self.initial_momentum_y = Some(diagnostics.momentum_y);
            self.initial_energy = Some(diagnostics.kinetic_energy);

            self.headers_written = true;
        }

        // Write data rows
        if let Err(e) = self.write_metrics(diagnostics) {
            warn!("Failed to write metrics: {}", e);
        }

        if let Err(e) = self.write_errors(diagnostics) {
            warn!("Failed to write errors: {}", e);
        }

        self.time_since_write = 0.0;
    }

    fn write_headers(&mut self) -> io::Result<()> {
        // Metrics file header
        if let Some(ref mut file) = self.metrics_file {
            writeln!(
                file,
                "frame,sim_time,mass_o2,mass_n2,mass_co2,momentum_x,momentum_y,kinetic_energy,\
                 max_velocity,max_divergence,avg_divergence,pressure_iters,timestep,max_cfl"
            )?;
            file.flush()?;
        }

        // Errors file header
        if let Some(ref mut file) = self.errors_file {
            writeln!(
                file,
                "frame,sim_time,mass_o2_error,mass_n2_error,mass_co2_error,\
                 momentum_x_error,momentum_y_error,energy_error"
            )?;
            file.flush()?;
        }

        Ok(())
    }

    fn write_metrics(&mut self, diagnostics: &SimulationDiagnostics) -> io::Result<()> {
        if let Some(ref mut file) = self.metrics_file {
            writeln!(
                file,
                "{},{:.3},{:.6},{:.6},{:.6},{:.3},{:.3},{:.3},{:.3},{:.6},{:.6},{},{:.6},{:.3}",
                diagnostics.frame,
                diagnostics.sim_time,
                diagnostics.mass_o2,
                diagnostics.mass_n2,
                diagnostics.mass_co2,
                diagnostics.momentum_x,
                diagnostics.momentum_y,
                diagnostics.kinetic_energy,
                diagnostics.max_velocity,
                diagnostics.max_divergence,
                diagnostics.avg_divergence,
                diagnostics.pressure_iterations,
                diagnostics.timestep,
                diagnostics.max_cfl,
            )?;
            file.flush()?;
        }
        Ok(())
    }

    fn write_errors(&mut self, diagnostics: &SimulationDiagnostics) -> io::Result<()> {
        if let Some(ref mut file) = self.errors_file {
            // Compute errors relative to initial values (or 0 if no initial value captured yet)
            let mass_o2_error = match self.initial_mass_o2 {
                Some(initial) => {
                    if initial > f32::EPSILON {
                        (diagnostics.mass_o2 - initial) / initial * 100.0
                    } else {
                        0.0
                    }
                }
                None => 0.0,
            };

            let mass_n2_error = match self.initial_mass_n2 {
                Some(initial) => {
                    if initial > f32::EPSILON {
                        (diagnostics.mass_n2 - initial) / initial * 100.0
                    } else {
                        0.0
                    }
                }
                None => 0.0,
            };

            let mass_co2_error = match self.initial_mass_co2 {
                Some(initial) => {
                    if initial > f32::EPSILON {
                        (diagnostics.mass_co2 - initial) / initial * 100.0
                    } else {
                        0.0
                    }
                }
                None => 0.0,
            };

            let momentum_x_error = match self.initial_momentum_x {
                Some(initial) => {
                    if initial.abs() > f32::EPSILON {
                        (diagnostics.momentum_x - initial) / initial * 100.0
                    } else {
                        diagnostics.momentum_x // Absolute error if started from zero
                    }
                }
                None => 0.0,
            };

            let momentum_y_error = match self.initial_momentum_y {
                Some(initial) => {
                    if initial.abs() > f32::EPSILON {
                        (diagnostics.momentum_y - initial) / initial * 100.0
                    } else {
                        diagnostics.momentum_y
                    }
                }
                None => 0.0,
            };

            let energy_error = match self.initial_energy {
                Some(initial) => {
                    if initial > f32::EPSILON {
                        (diagnostics.kinetic_energy - initial) / initial * 100.0
                    } else {
                        0.0
                    }
                }
                None => 0.0,
            };

            writeln!(
                file,
                "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}",
                diagnostics.frame,
                diagnostics.sim_time,
                mass_o2_error,
                mass_n2_error,
                mass_co2_error,
                momentum_x_error,
                momentum_y_error,
                energy_error,
            )?;
            file.flush()?;
        }
        Ok(())
    }

    /// Flush and close files (call this on shutdown if needed).
    pub fn finalize(&mut self) {
        if let Some(ref mut file) = self.metrics_file {
            let _ = file.flush();
        }
        if let Some(ref mut file) = self.errors_file {
            let _ = file.flush();
        }
        self.metrics_file = None;
        self.errors_file = None;
    }
}

impl Drop for CsvExporter {
    fn drop(&mut self) {
        self.finalize();
    }
}

/// Estimate file size for a simulation run.
/// Returns approximate size in MB for metrics and errors files combined.
#[allow(dead_code)]
pub fn estimate_file_size(duration_seconds: f32, write_interval: f32) -> f32 {
    let num_entries = (duration_seconds / write_interval).ceil();
    let bytes_per_metrics_row = 100.0; // Approximate bytes per CSV row
    let bytes_per_errors_row = 80.0;
    let total_bytes = num_entries * (bytes_per_metrics_row + bytes_per_errors_row);
    total_bytes / 1_000_000.0 // Convert to MB
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_size_estimate() {
        // 1 hour at 5 second intervals = 720 entries
        let size_mb = estimate_file_size(3600.0, 5.0);
        // Should be around 0.13 MB (very small)
        assert!(size_mb < 1.0, "1 hour should be less than 1 MB");

        // 24 hours at 5 second intervals = 17,280 entries
        let size_mb = estimate_file_size(86400.0, 5.0);
        // Should be around 3.1 MB (still very reasonable)
        assert!(size_mb < 10.0, "24 hours should be less than 10 MB");
    }
}
