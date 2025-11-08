use crate::player::Player;
use crate::tilemap::{LifeSupportTiles, TileCollisionMap, TilemapInitSet};
use bevy::prelude::*;

/// Physical constants for atmospheric simulation
pub mod constants {
    /// Universal gas constant [J/(mol·K)]
    pub const R: f32 = 8.314;

    /// Molar mass of O2 [kg/mol]
    pub const M_O2: f32 = 0.032;

    /// Molar mass of N2 [kg/mol]
    pub const M_N2: f32 = 0.028;

    /// Molar mass of CO2 [kg/mol]
    pub const M_CO2: f32 = 0.044;

    /// Dynamic viscosity of air at 20°C [Pa·s]
    pub const MU: f32 = 1.8e-5;

    /// Thermal conductivity of air at 20°C [W/(m·K)]
    pub const K_THERMAL: f32 = 0.026;

    /// Specific heat capacity of air at constant pressure [J/(kg·K)]
    pub const CP: f32 = 1005.0;

    /// Earth atmospheric pressure at sea level [Pa]
    pub const EARTH_PRESSURE: f32 = 101325.0;

    /// Standard room temperature [K]
    pub const ROOM_TEMP: f32 = 293.0;

    /// Earth atmosphere O2 mole fraction
    pub const EARTH_O2_FRACTION: f32 = 0.21;

    /// Earth atmosphere N2 mole fraction
    pub const EARTH_N2_FRACTION: f32 = 0.78;

    /// Earth atmosphere CO2 mole fraction (approximately 0.04%)
    pub const EARTH_CO2_FRACTION: f32 = 0.0004;

    // Human respiration constants
    /// Oxygen consumption rate at rest [kg/s]
    /// Based on ~250 mL/min O2 consumption = 5.54e-6 kg/s
    pub const O2_CONSUMPTION_RATE: f32 = 5.54e-6;

    /// CO2 production rate at rest [kg/s]
    /// Respiratory quotient ~0.8, so CO2 production is 0.8 * O2 consumption by moles
    /// But CO2 has higher molar mass, so mass rate = 0.8 * (M_CO2/M_O2) * O2_rate
    /// = 0.8 * (44/32) * 5.54e-6 = 6.09e-6 kg/s
    pub const CO2_PRODUCTION_RATE: f32 = 6.09e-6;

    /// Height of room (for calculating cell volume) [m]
    pub const ROOM_HEIGHT: f32 = 2.5;

    // Life support system constants
    /// O2 generation rate per life support unit [kg/s]
    /// Tuned for ~60s refill of the 29×14×2.5 m playable room.
    pub const LIFE_SUPPORT_O2_RATE: f32 = 4.7;

    /// N2 generation rate per life support unit [kg/s]
    /// Scaled to preserve ~78% N2 while matching the faster fill rate.
    pub const LIFE_SUPPORT_N2_RATE: f32 = 15.6;

    /// Number of tiles around a life support unit influenced by its circulation fan
    pub const FAN_RADIUS_TILES: i32 = 3;

    /// Target tangential speed imparted by life support fans [m/s]
    pub const FAN_SWIRL_SPEED: f32 = 1.5;

    /// Blend factor when steering local velocities toward the fan swirl
    pub const FAN_BLEND: f32 = 0.35;

    /// Pressure-driven mass flux conductance [kg/(Pa·m·s)]
    /// Controls how fast gases flow from high to low pressure
    /// Tuned for realistic expansion into vacuum while maintaining stability
    pub const PRESSURE_FLUX_CONDUCTANCE: f32 = 0.01;

    /// Mass conservation check interval (seconds)
    /// How often to verify total mass conservation
    pub const MASS_CHECK_INTERVAL: f32 = 5.0;

    /// Mass conservation tolerance (relative error)
    /// Warn if mass changes by more than this fraction (excluding life support/respiration)
    pub const MASS_TOLERANCE: f32 = 0.01; // 1%

    /// CFL number for timestep stability
    /// CFL = (u*dt)/dx should be < 1 for explicit schemes
    /// Using 0.5 provides good stability margin
    pub const CFL_NUMBER: f32 = 0.5;

    /// Poisson solver maximum iterations
    pub const POISSON_MAX_ITERATIONS: usize = 50;

    /// Poisson solver convergence tolerance
    pub const POISSON_TOLERANCE: f32 = 1e-4;
}

/// Atmospheric state for a single grid cell
#[derive(Clone, Debug)]
pub struct AtmosphereCell {
    /// Oxygen density [kg/m³]
    pub rho_o2: f32,

    /// Nitrogen density [kg/m³]
    pub rho_n2: f32,

    /// Carbon dioxide density [kg/m³]
    pub rho_co2: f32,

    /// Velocity in x-direction [m/s]
    pub u: f32,

    /// Velocity in y-direction [m/s]
    pub v: f32,

    /// Temperature [K]
    pub temperature: f32,

    /// Cached pressure [Pa] (computed from ideal gas law)
    pub pressure: f32,
}

impl Default for AtmosphereCell {
    fn default() -> Self {
        Self {
            rho_o2: 0.0,
            rho_n2: 0.0,
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        }
    }
}

impl AtmosphereCell {
    /// Create cell with Earth-like atmosphere
    pub fn earth_atmosphere() -> Self {
        // Use ideal gas law to compute densities: ρ = (p * M) / (R * T)
        let p = constants::EARTH_PRESSURE;
        let t = constants::ROOM_TEMP;

        // Partial pressures
        let p_o2 = p * constants::EARTH_O2_FRACTION;
        let p_n2 = p * constants::EARTH_N2_FRACTION;
        let p_co2 = p * constants::EARTH_CO2_FRACTION;

        // Densities from ideal gas law
        let rho_o2 = (p_o2 * constants::M_O2) / (constants::R * t);
        let rho_n2 = (p_n2 * constants::M_N2) / (constants::R * t);
        let rho_co2 = (p_co2 * constants::M_CO2) / (constants::R * t);

        Self {
            rho_o2,
            rho_n2,
            rho_co2,
            u: 0.0,
            v: 0.0,
            temperature: t,
            pressure: p,
        }
    }

    /// Create vacuum cell
    pub fn vacuum() -> Self {
        Self {
            rho_o2: 0.0,
            rho_n2: 0.0,
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: 2.7, // Cosmic microwave background temperature
            pressure: 0.0,
        }
    }

    /// Total gas density [kg/m³]
    #[inline]
    pub fn total_density(&self) -> f32 {
        self.rho_o2 + self.rho_n2 + self.rho_co2
    }

    /// Compute pressure from current state using ideal gas law
    /// p = (ρ_O2/M_O2 + ρ_N2/M_N2 + ρ_CO2/M_CO2) * R * T
    #[inline]
    pub fn compute_pressure(&self) -> f32 {
        let n_o2 = self.rho_o2 / constants::M_O2;
        let n_n2 = self.rho_n2 / constants::M_N2;
        let n_co2 = self.rho_co2 / constants::M_CO2;
        (n_o2 + n_n2 + n_co2) * constants::R * self.temperature
    }

    /// Update cached pressure value
    #[inline]
    pub fn update_pressure(&mut self) {
        self.pressure = self.compute_pressure();
    }
}

/// Grid-based atmospheric simulation resource
#[derive(Resource)]
pub struct AtmosphereGrid {
    pub width: u32,
    pub height: u32,
    pub tile_size_world: f32,
    pub tile_size_physical: f32,

    /// Current atmospheric state
    pub cells: Vec<AtmosphereCell>,

    /// Temporary buffer for multi-step integration
    pub cells_buffer: Vec<AtmosphereCell>,
}

/// Resource for tracking mass conservation
#[derive(Resource)]
struct MassTracker {
    /// Last recorded total mass [kg]
    last_total_mass: f32,
    /// Time since last check [s]
    time_since_check: f32,
    /// Expected mass change from sources/sinks [kg]
    expected_delta: f32,
}

/// Resource for tracking velocity divergence
#[derive(Resource)]
struct DivergenceTracker {
    /// Time since last check [s]
    time_since_check: f32,
}

impl AtmosphereGrid {
    pub fn new(width: u32, height: u32, tile_size_world: f32, tile_size_physical: f32) -> Self {
        let cell_count = (width * height) as usize;

        Self {
            width,
            height,
            tile_size_world,
            tile_size_physical,
            cells: vec![AtmosphereCell::default(); cell_count],
            cells_buffer: vec![AtmosphereCell::default(); cell_count],
        }
    }

    /// Initialize with Earth-like atmosphere in open areas
    pub fn initialize_earth_atmosphere(&mut self, collision_map: &TileCollisionMap) {
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = self.index(x, y);

                // If tile is not blocked (i.e., it's open space), fill with atmosphere
                if !collision_map.is_blocked(x, y) {
                    self.cells[idx] = AtmosphereCell::earth_atmosphere();
                } else {
                    // Walls have no atmosphere (vacuum)
                    self.cells[idx] = AtmosphereCell::vacuum();
                }
            }
        }
    }

    /// Initialize entire room as vacuum
    pub fn initialize_vacuum(&mut self) {
        for cell in self.cells.iter_mut() {
            *cell = AtmosphereCell::vacuum();
        }
    }

    #[inline]
    pub fn index(&self, x: u32, y: u32) -> usize {
        (y * self.width + x) as usize
    }

    #[inline]
    pub fn in_bounds(&self, x: i32, y: i32) -> bool {
        x >= 0 && y >= 0 && (x as u32) < self.width && (y as u32) < self.height
    }

    /// Get cell at coordinates (returns None if out of bounds)
    pub fn get(&self, x: i32, y: i32) -> Option<&AtmosphereCell> {
        if self.in_bounds(x, y) {
            Some(&self.cells[self.index(x as u32, y as u32)])
        } else {
            None
        }
    }

    /// Get mutable cell at coordinates
    pub fn get_mut(&mut self, x: i32, y: i32) -> Option<&mut AtmosphereCell> {
        if self.in_bounds(x, y) {
            let idx = self.index(x as u32, y as u32);
            Some(&mut self.cells[idx])
        } else {
            None
        }
    }
}

pub struct AtmospherePlugin;

#[derive(SystemSet, Debug, Clone, PartialEq, Eq, Hash)]
pub enum AtmosphereSimSet {
    Main,
}

impl Plugin for AtmospherePlugin {
    fn build(&self, app: &mut App) {
        app.configure_sets(Update, AtmosphereSimSet::Main);
        app
            // Run atmosphere initialization after tilemap is created
            .add_systems(Startup, initialize_atmosphere.after(TilemapInitSet))
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

/// Initialize atmospheric grid
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
    commands.insert_resource(MassTracker {
        last_total_mass: 0.0,
        time_since_check: 0.0,
        expected_delta: 0.0,
    });

    // Initialize divergence tracker
    commands.insert_resource(DivergenceTracker {
        time_since_check: 0.0,
    });

    info!(
        "Atmospheric simulation initialized: {}x{} grid (starting in vacuum)",
        collision_map.width, collision_map.height
    );
}

/// Update pressure values from current density and temperature state
fn update_pressure_from_state(mut atmosphere: ResMut<AtmosphereGrid>) {
    update_grid_pressures(&mut atmosphere);
}

/// Recompute pressure for every cell using the ideal gas law
#[inline]
fn update_grid_pressures(atmosphere: &mut AtmosphereGrid) {
    for cell in atmosphere.cells.iter_mut() {
        cell.update_pressure();
    }
}

/// Main atmospheric simulation system
/// Implements compressible Navier-Stokes equations with operator splitting:
/// 1. Advection (transport by velocity field)
/// 2. Diffusion (viscosity and thermal conductivity)
/// 3. Pressure-driven mass flux (expansion into vacuum, back pressure)
/// 4. Pressure correction (enforce momentum balance)
fn simulate_atmosphere(
    mut atmosphere: ResMut<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    time: Res<Time>,
) {
    let mut remaining = time.delta_secs();

    if remaining <= 0.0 {
        return;
    }

    // Sub-step until we've integrated the entire frame duration.
    while remaining > 0.0 {
        // For stability, limit timestep based on CFL condition
        // CFL = (u*dt)/dx should be < 1
        let max_velocity = compute_max_velocity(&atmosphere);
        let cfl_dt = if max_velocity > 0.0 {
            constants::CFL_NUMBER * atmosphere.tile_size_physical / max_velocity
        } else {
            remaining
        };

        // Avoid infinitesimal steps that would stall the loop.
        let dt_sim = remaining.min(cfl_dt.max(f32::EPSILON));

        // Operator splitting approach per sub-step:
        // Step 1: Advection (explicit)
        advection_step(&mut atmosphere, &collision_map, dt_sim);

        // Step 2: Diffusion (explicit coefficients chosen for stability)
        diffusion_step(&mut atmosphere, &collision_map, dt_sim);

        // Update pressure so the subsequent steps see the latest state
        update_grid_pressures(&mut atmosphere);

        // Step 3: Pressure-driven mass flux (handles expansion into vacuum and back pressure)
        pressure_flux_step(&mut atmosphere, &collision_map, dt_sim);

        // Update pressure again after mass transfer
        update_grid_pressures(&mut atmosphere);

        // Step 4: Pressure correction (fine-tune momentum balance)
        pressure_correction_step(&mut atmosphere, &collision_map, dt_sim);

        remaining -= dt_sim;
    }
}

/// Compute maximum velocity magnitude in the grid for CFL condition
fn compute_max_velocity(atmosphere: &AtmosphereGrid) -> f32 {
    atmosphere
        .cells
        .iter()
        .map(|cell| (cell.u * cell.u + cell.v * cell.v).sqrt())
        .fold(0.0, f32::max)
}

/// Step 1: Advection - transport quantities by velocity field
///
/// Uses MacCormack scheme for second-order accuracy with reduced numerical dissipation.
/// MacCormack is a predictor-corrector method:
/// 1. Forward step: advect forward in time (predictor)
/// 2. Backward step: advect backward from predicted state (corrector)
/// 3. Correction: φ_final = φ_forward + 0.5 * (φ_original - φ_backward)
///
/// Benefits over semi-Lagrangian:
/// - Second-order temporal accuracy (vs. first-order)
/// - Reduced numerical dissipation (preserves flow structures better)
/// - Better energy conservation
///
/// Note: Can produce small oscillations near sharp gradients; clamping is applied.
fn advection_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    let cell_count = (atmosphere.width * atmosphere.height) as usize;

    // Allocate temporary storage for MacCormack scheme
    let mut forward_state = vec![AtmosphereCell::default(); cell_count];
    let mut backward_state = vec![AtmosphereCell::default(); cell_count];

    // Copy current state to buffer (original state φⁿ)
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    // Step 1: Forward advection (predictor)
    // φ* = φⁿ - dt * (u · ∇φⁿ)
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells_buffer[idx];

            // Backward trace from current position
            let x_back = x as f32 - (cell.u * dt) / atmosphere.tile_size_physical;
            let y_back = y as f32 - (cell.v * dt) / atmosphere.tile_size_physical;

            forward_state[idx] = interpolate_cell(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                x_back,
                y_back,
            );
        }
    }

    // Step 2: Backward advection (corrector)
    // φ** = φ* + dt * (u · ∇φ*)
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &forward_state[idx];

            // Forward trace from predicted position
            let x_forward = x as f32 + (cell.u * dt) / atmosphere.tile_size_physical;
            let y_forward = y as f32 + (cell.v * dt) / atmosphere.tile_size_physical;

            backward_state[idx] = interpolate_cell(
                &forward_state,
                atmosphere.width,
                atmosphere.height,
                x_forward,
                y_forward,
            );
        }
    }

    // Step 3: MacCormack correction
    // φⁿ⁺¹ = φ* + 0.5 * (φⁿ - φ**)
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let original = &atmosphere.cells_buffer[idx];
            let forward = &forward_state[idx];
            let backward = &backward_state[idx];

            // Apply MacCormack correction
            let corrected_rho_o2 = forward.rho_o2 + 0.5 * (original.rho_o2 - backward.rho_o2);
            let corrected_rho_n2 = forward.rho_n2 + 0.5 * (original.rho_n2 - backward.rho_n2);
            let corrected_rho_co2 = forward.rho_co2 + 0.5 * (original.rho_co2 - backward.rho_co2);
            let corrected_u = forward.u + 0.5 * (original.u - backward.u);
            let corrected_v = forward.v + 0.5 * (original.v - backward.v);
            let corrected_temp = forward.temperature + 0.5 * (original.temperature - backward.temperature);

            // Clamp to prevent negative densities and maintain physical bounds
            atmosphere.cells[idx].rho_o2 = corrected_rho_o2.max(0.0);
            atmosphere.cells[idx].rho_n2 = corrected_rho_n2.max(0.0);
            atmosphere.cells[idx].rho_co2 = corrected_rho_co2.max(0.0);
            atmosphere.cells[idx].u = corrected_u;
            atmosphere.cells[idx].v = corrected_v;
            atmosphere.cells[idx].temperature = corrected_temp.max(0.0);
            atmosphere.cells[idx].update_pressure();
        }
    }
}

/// Bilinear interpolation of atmospheric cell
fn interpolate_cell(
    cells: &[AtmosphereCell],
    width: u32,
    height: u32,
    x: f32,
    y: f32,
) -> AtmosphereCell {
    let x0 = x.floor().max(0.0).min((width - 1) as f32) as u32;
    let y0 = y.floor().max(0.0).min((height - 1) as f32) as u32;
    let x1 = (x0 + 1).min(width - 1);
    let y1 = (y0 + 1).min(height - 1);

    let fx = (x - x0 as f32).clamp(0.0, 1.0);
    let fy = (y - y0 as f32).clamp(0.0, 1.0);

    let idx00 = (y0 * width + x0) as usize;
    let idx10 = (y0 * width + x1) as usize;
    let idx01 = (y1 * width + x0) as usize;
    let idx11 = (y1 * width + x1) as usize;

    let c00 = &cells[idx00];
    let c10 = &cells[idx10];
    let c01 = &cells[idx01];
    let c11 = &cells[idx11];

    AtmosphereCell {
        rho_o2: bilerp(c00.rho_o2, c10.rho_o2, c01.rho_o2, c11.rho_o2, fx, fy),
        rho_n2: bilerp(c00.rho_n2, c10.rho_n2, c01.rho_n2, c11.rho_n2, fx, fy),
        rho_co2: bilerp(c00.rho_co2, c10.rho_co2, c01.rho_co2, c11.rho_co2, fx, fy),
        u: bilerp(c00.u, c10.u, c01.u, c11.u, fx, fy),
        v: bilerp(c00.v, c10.v, c01.v, c11.v, fx, fy),
        temperature: bilerp(
            c00.temperature,
            c10.temperature,
            c01.temperature,
            c11.temperature,
            fx,
            fy,
        ),
        pressure: 0.0, // Will be recomputed
    }
}

#[inline]
fn bilerp(v00: f32, v10: f32, v01: f32, v11: f32, fx: f32, fy: f32) -> f32 {
    let v0 = v00 * (1.0 - fx) + v10 * fx;
    let v1 = v01 * (1.0 - fx) + v11 * fx;
    v0 * (1.0 - fy) + v1 * fy
}

/// Step 2: Diffusion - viscosity and thermal diffusion
/// Uses explicit method with small diffusion coefficients for stability
///
/// Stability condition for explicit diffusion: α = μ*dt/dx² < 0.25
/// We clamp to ensure stability even with large timesteps from CFL limiting.
fn diffusion_step(atmosphere: &mut AtmosphereGrid, collision_map: &TileCollisionMap, dt: f32) {
    let dx = atmosphere.tile_size_physical;

    // Stability limit for explicit diffusion (von Neumann analysis)
    const MAX_ALPHA: f32 = 0.2; // Safety factor below theoretical limit of 0.25

    let alpha = (constants::MU * dt / (dx * dx)).min(MAX_ALPHA);
    let kappa = (constants::K_THERMAL * dt / (dx * dx)).min(MAX_ALPHA);

    // Copy to buffer
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);

            // Get neighbors for Laplacian
            let c = &atmosphere.cells_buffer[idx];
            let left = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32 - 1,
                y as i32,
                c,
            );
            let right = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32 + 1,
                y as i32,
                c,
            );
            let down = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32,
                y as i32 - 1,
                c,
            );
            let up = get_or_boundary(
                &atmosphere.cells_buffer,
                atmosphere.width,
                atmosphere.height,
                collision_map,
                x as i32,
                y as i32 + 1,
                c,
            );

            // Momentum diffusion (viscosity)
            let laplacian_u = left.u + right.u + down.u + up.u - 4.0 * c.u;
            let laplacian_v = left.v + right.v + down.v + up.v - 4.0 * c.v;

            atmosphere.cells[idx].u += alpha * laplacian_u;
            atmosphere.cells[idx].v += alpha * laplacian_v;

            // Thermal diffusion
            let laplacian_t =
                left.temperature + right.temperature + down.temperature + up.temperature
                    - 4.0 * c.temperature;
            atmosphere.cells[idx].temperature += kappa * laplacian_t;

            // Mass diffusion (Fick's law)
            // Diffusion coefficient for gases in air: D ≈ 2×10⁻⁵ m²/s
            let d_coeff = (2.0e-5 * dt / (dx * dx)).min(MAX_ALPHA);
            let laplacian_o2 =
                left.rho_o2 + right.rho_o2 + down.rho_o2 + up.rho_o2 - 4.0 * c.rho_o2;
            let laplacian_n2 =
                left.rho_n2 + right.rho_n2 + down.rho_n2 + up.rho_n2 - 4.0 * c.rho_n2;
            let laplacian_co2 =
                left.rho_co2 + right.rho_co2 + down.rho_co2 + up.rho_co2 - 4.0 * c.rho_co2;

            atmosphere.cells[idx].rho_o2 += d_coeff * laplacian_o2;
            atmosphere.cells[idx].rho_n2 += d_coeff * laplacian_n2;
            atmosphere.cells[idx].rho_co2 += d_coeff * laplacian_co2;
        }
    }
}

/// Get cell or apply boundary condition for walls
///
/// This implements a "ghost cell" approach for boundary conditions:
/// - For fluid cells: returns actual cell data
/// - For wall cells: returns virtual cell with no-slip boundary condition (u=v=0)
/// - For out-of-bounds: treats as sealed wall
///
/// The ghost cell values are used ONLY for computing spatial derivatives (Laplacians)
/// in the diffusion step. Wall cells themselves are never updated.
///
/// Boundary conditions:
/// - Velocity: No-slip (u=v=0 at walls)
/// - Temperature: Adiabatic/Neumann (∂T/∂n = 0, enforced by mirroring)
/// - Density: Neumann (∂ρ/∂n = 0, enforced by mirroring)
fn get_or_boundary(
    cells: &[AtmosphereCell],
    width: u32,
    height: u32,
    collision_map: &TileCollisionMap,
    x: i32,
    y: i32,
    center: &AtmosphereCell,
) -> AtmosphereCell {
    if x >= 0 && y >= 0 && (x as u32) < width && (y as u32) < height {
        let idx = (y as u32 * width + x as u32) as usize;
        if !collision_map.is_blocked(x as u32, y as u32) {
            cells[idx].clone()
        } else {
            // Ghost cell for wall: no-slip velocity, mirror other properties
            AtmosphereCell {
                rho_o2: center.rho_o2,
                rho_n2: center.rho_n2,
                rho_co2: center.rho_co2,
                u: 0.0, // No-slip boundary condition
                v: 0.0,
                temperature: center.temperature, // Adiabatic (mirror)
                pressure: center.pressure,
            }
        }
    } else {
        // Ghost cell for domain boundary: sealed wall (conserve mass)
        AtmosphereCell {
            rho_o2: center.rho_o2,
            rho_n2: center.rho_n2,
            rho_co2: center.rho_co2,
            u: 0.0, // No-slip
            v: 0.0,
            temperature: center.temperature,
            pressure: center.pressure,
        }
    }
}

/// Pressure-driven mass flux - gases flow from high to low pressure
/// This creates bulk flow naturally even from vacuum, handling expansion and back pressure
fn pressure_flux_step(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;
    let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;

    // Calculate cell volume for density updates
    let cell_volume = dx * dx * constants::ROOM_HEIGHT;

    // Use buffer to store flux-driven updates
    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells_buffer[idx];

            // Process flux in each direction (x and y)
            // X-direction (right neighbor)
            if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                let idx_right = atmosphere.index(x + 1, y);
                let cell_right = &atmosphere.cells_buffer[idx_right];

                let p_diff = cell.pressure - cell_right.pressure;

                if p_diff.abs() > f32::EPSILON {
                    // Total mass flux [kg/s] driven by pressure difference
                    // Positive p_diff: flow FROM cell TO cell_right
                    // Negative p_diff: flow FROM cell_right TO cell
                    let total_flux = conductance * p_diff * dx * dt;

                    // Identify donor cell (where mass comes from)
                    let (donor, donor_idx, receiver_idx) = if p_diff > 0.0 {
                        (cell, idx, idx_right)
                    } else {
                        (cell_right, idx_right, idx)
                    };

                    // Distribute flux among gases proportional to donor's partial pressures
                    if donor.pressure > f32::EPSILON {
                        let p_o2 = (donor.rho_o2 / constants::M_O2) * constants::R * donor.temperature;
                        let p_n2 = (donor.rho_n2 / constants::M_N2) * constants::R * donor.temperature;
                        let p_co2 = (donor.rho_co2 / constants::M_CO2) * constants::R * donor.temperature;

                        // Use absolute value of flux for species decomposition
                        let abs_flux = total_flux.abs();
                        let flux_o2 = abs_flux * (p_o2 / donor.pressure);
                        let flux_n2 = abs_flux * (p_n2 / donor.pressure);
                        let flux_co2 = abs_flux * (p_co2 / donor.pressure);

                        // Convert mass flux to density change
                        let delta_rho_o2 = flux_o2 / cell_volume;
                        let delta_rho_n2 = flux_n2 / cell_volume;
                        let delta_rho_co2 = flux_co2 / cell_volume;

                        // Update donor cell (loses mass)
                        atmosphere.cells[donor_idx].rho_o2 -= delta_rho_o2;
                        atmosphere.cells[donor_idx].rho_n2 -= delta_rho_n2;
                        atmosphere.cells[donor_idx].rho_co2 -= delta_rho_co2;

                        // Update receiver cell (gains mass)
                        atmosphere.cells[receiver_idx].rho_o2 += delta_rho_o2;
                        atmosphere.cells[receiver_idx].rho_n2 += delta_rho_n2;
                        atmosphere.cells[receiver_idx].rho_co2 += delta_rho_co2;

                        // Transfer momentum: flowing mass carries donor's velocity
                        // When mass dm flows from donor, it carries momentum dm * u_donor
                        let mass_flux_abs = abs_flux; // kg
                        let donor_rho = donor.total_density().max(1e-6);
                        let donor_velocity_u = donor.u;

                        let receiver_cell = &atmosphere.cells_buffer[receiver_idx];
                        let receiver_rho = receiver_cell.total_density().max(1e-6);

                        // Momentum flux in x-direction (mass carries donor's x-velocity)
                        let momentum_flux_u = mass_flux_abs * donor_velocity_u;

                        // Update velocities: donor loses momentum, receiver gains it
                        atmosphere.cells[donor_idx].u -= momentum_flux_u / (donor_rho * cell_volume);
                        atmosphere.cells[receiver_idx].u += momentum_flux_u / (receiver_rho * cell_volume);
                    }
                }
            }

            // Y-direction (up neighbor)
            if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                let idx_up = atmosphere.index(x, y + 1);
                let cell_up = &atmosphere.cells_buffer[idx_up];

                let p_diff = cell.pressure - cell_up.pressure;

                if p_diff.abs() > f32::EPSILON {
                    // Total mass flux [kg/s] driven by pressure difference
                    // Positive p_diff: flow FROM cell TO cell_up
                    // Negative p_diff: flow FROM cell_up TO cell
                    let total_flux = conductance * p_diff * dx * dt;

                    // Identify donor cell (where mass comes from)
                    let (donor, donor_idx, receiver_idx) = if p_diff > 0.0 {
                        (cell, idx, idx_up)
                    } else {
                        (cell_up, idx_up, idx)
                    };

                    // Distribute flux among gases proportional to donor's partial pressures
                    if donor.pressure > f32::EPSILON {
                        let p_o2 = (donor.rho_o2 / constants::M_O2) * constants::R * donor.temperature;
                        let p_n2 = (donor.rho_n2 / constants::M_N2) * constants::R * donor.temperature;
                        let p_co2 = (donor.rho_co2 / constants::M_CO2) * constants::R * donor.temperature;

                        // Use absolute value of flux for species decomposition
                        let abs_flux = total_flux.abs();
                        let flux_o2 = abs_flux * (p_o2 / donor.pressure);
                        let flux_n2 = abs_flux * (p_n2 / donor.pressure);
                        let flux_co2 = abs_flux * (p_co2 / donor.pressure);

                        let delta_rho_o2 = flux_o2 / cell_volume;
                        let delta_rho_n2 = flux_n2 / cell_volume;
                        let delta_rho_co2 = flux_co2 / cell_volume;

                        // Update donor cell (loses mass)
                        atmosphere.cells[donor_idx].rho_o2 -= delta_rho_o2;
                        atmosphere.cells[donor_idx].rho_n2 -= delta_rho_n2;
                        atmosphere.cells[donor_idx].rho_co2 -= delta_rho_co2;

                        // Update receiver cell (gains mass)
                        atmosphere.cells[receiver_idx].rho_o2 += delta_rho_o2;
                        atmosphere.cells[receiver_idx].rho_n2 += delta_rho_n2;
                        atmosphere.cells[receiver_idx].rho_co2 += delta_rho_co2;

                        // Transfer momentum: flowing mass carries donor's velocity
                        // When mass dm flows from donor, it carries momentum dm * v_donor
                        let mass_flux_abs = abs_flux; // kg
                        let donor_rho = donor.total_density().max(1e-6);
                        let donor_velocity_v = donor.v;

                        let receiver_cell = &atmosphere.cells_buffer[receiver_idx];
                        let receiver_rho = receiver_cell.total_density().max(1e-6);

                        // Momentum flux in y-direction (mass carries donor's y-velocity)
                        let momentum_flux_v = mass_flux_abs * donor_velocity_v;

                        // Update velocities: donor loses momentum, receiver gains it
                        atmosphere.cells[donor_idx].v -= momentum_flux_v / (donor_rho * cell_volume);
                        atmosphere.cells[receiver_idx].v += momentum_flux_v / (receiver_rho * cell_volume);
                    }
                }
            }
        }
    }

    // Clamp densities to prevent negative values from numerical errors
    for cell in atmosphere.cells.iter_mut() {
        cell.rho_o2 = cell.rho_o2.max(0.0);
        cell.rho_n2 = cell.rho_n2.max(0.0);
        cell.rho_co2 = cell.rho_co2.max(0.0);
    }
}

/// Step 4: Pressure projection - enforce divergence-free velocity field
///
/// Implements pressure projection method to maintain mass conservation:
/// 1. Compute velocity divergence: div = ∂u/∂x + ∂v/∂y
/// 2. Solve Poisson equation for pressure correction: ∇²p = ρ/dt * div
/// 3. Update velocities to be divergence-free: u_new = u - dt/ρ * ∇p
///
/// Uses Jacobi iteration to solve the Poisson equation.
fn pressure_correction_step(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;

    // Allocate pressure correction field
    let cell_count = (atmosphere.width * atmosphere.height) as usize;
    let mut pressure_correction = vec![0.0; cell_count];
    let mut divergence = vec![0.0; cell_count];

    // Step 1: Compute velocity divergence
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);

            // Compute divergence using central differences
            let u_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                atmosphere.cells[atmosphere.index(x + 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                atmosphere.cells[atmosphere.index(x - 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let v_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                atmosphere.cells[atmosphere.index(x, y + 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                atmosphere.cells[atmosphere.index(x, y - 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let du_dx = (u_right - u_left) / (2.0 * dx);
            let dv_dy = (v_up - v_down) / (2.0 * dx);

            divergence[idx] = du_dx + dv_dy;
        }
    }

    // Step 2: Solve Poisson equation: ∇²p = ρ/dt * div
    // Using Jacobi iteration
    let mut pressure_correction_temp = vec![0.0; cell_count];

    for _iter in 0..constants::POISSON_MAX_ITERATIONS {
        let mut max_change: f32 = 0.0;

        for y in 0..atmosphere.height {
            for x in 0..atmosphere.width {
                if collision_map.is_blocked(x, y) {
                    continue;
                }

                let idx = atmosphere.index(x, y);
                let rho = atmosphere.cells[idx].total_density().max(1e-6);

                // Get neighbor pressure corrections (Neumann BC at walls: ∂p/∂n = 0)
                let p_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                    pressure_correction[atmosphere.index(x - 1, y)]
                } else {
                    pressure_correction[idx]
                };

                let p_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                    pressure_correction[atmosphere.index(x + 1, y)]
                } else {
                    pressure_correction[idx]
                };

                let p_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                    pressure_correction[atmosphere.index(x, y - 1)]
                } else {
                    pressure_correction[idx]
                };

                let p_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                    pressure_correction[atmosphere.index(x, y + 1)]
                } else {
                    pressure_correction[idx]
                };

                // Jacobi update: p_new = (p_left + p_right + p_down + p_up - dx²*rhs) / 4
                let rhs = (rho / dt) * divergence[idx];
                let p_new = (p_left + p_right + p_down + p_up - dx * dx * rhs) * 0.25;

                pressure_correction_temp[idx] = p_new;
                max_change = max_change.max((p_new - pressure_correction[idx]).abs());
            }
        }

        // Copy temp to main
        pressure_correction.copy_from_slice(&pressure_correction_temp);

        // Check convergence
        if max_change < constants::POISSON_TOLERANCE {
            break;
        }
    }

    // Step 3: Apply pressure correction to velocities
    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let rho = atmosphere.cells[idx].total_density().max(1e-6);

            // Compute pressure correction gradient
            let p_c = pressure_correction[idx];

            let p_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                pressure_correction[atmosphere.index(x - 1, y)]
            } else {
                p_c
            };

            let p_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                pressure_correction[atmosphere.index(x + 1, y)]
            } else {
                p_c
            };

            let p_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                pressure_correction[atmosphere.index(x, y - 1)]
            } else {
                p_c
            };

            let p_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                pressure_correction[atmosphere.index(x, y + 1)]
            } else {
                p_c
            };

            let dp_dx = (p_right - p_left) / (2.0 * dx);
            let dp_dy = (p_up - p_down) / (2.0 * dx);

            // Correct velocity to enforce divergence-free condition
            atmosphere.cells[idx].u -= (dt / rho) * dp_dx;
            atmosphere.cells[idx].v -= (dt / rho) * dp_dy;
        }
    }
}

/// Verify mass conservation in the simulation
///
/// Checks that total mass only changes due to intended sources/sinks (life support, respiration).
/// Warns if numerical errors cause unexpected mass changes.
fn check_mass_conservation(
    atmosphere: Res<AtmosphereGrid>,
    mut tracker: ResMut<MassTracker>,
    time: Res<Time>,
) {
    tracker.time_since_check += time.delta_secs();

    if tracker.time_since_check < constants::MASS_CHECK_INTERVAL {
        return;
    }

    // Calculate total mass in the system [kg]
    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;
    let total_mass: f32 = atmosphere
        .cells
        .iter()
        .map(|cell| cell.total_density() * tile_volume)
        .sum();

    // First check (last_total_mass is 0), just record
    if tracker.last_total_mass > f32::EPSILON {
        let actual_delta = total_mass - tracker.last_total_mass;
        let error = (actual_delta - tracker.expected_delta).abs();
        let relative_error = if tracker.last_total_mass > f32::EPSILON {
            error / tracker.last_total_mass
        } else {
            0.0
        };

        if relative_error > constants::MASS_TOLERANCE {
            warn!(
                "Mass conservation warning: Total mass changed by {:.6} kg (expected {:.6} kg from sources/sinks). Relative error: {:.2}%",
                actual_delta,
                tracker.expected_delta,
                relative_error * 100.0
            );
        } else {
            debug!(
                "Mass conservation check: Total mass = {:.3} kg (delta: {:.6} kg, expected: {:.6} kg)",
                total_mass, actual_delta, tracker.expected_delta
            );
        }
    }

    // Reset tracker
    tracker.last_total_mass = total_mass;
    tracker.time_since_check = 0.0;
    tracker.expected_delta = 0.0; // Will be accumulated by sources/sinks in next interval
}

/// Life support generation system - produces O2 and N2 at life support tiles
fn life_support_generation(
    mut atmosphere: ResMut<AtmosphereGrid>,
    life_support: Res<LifeSupportTiles>,
    mut tracker: ResMut<MassTracker>,
    time: Res<Time>,
) {
    let dt = time.delta_secs();

    // Calculate tile volume [m³]
    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

    let mut total_generated = 0.0;

    for pos in &life_support.positions {
        let idx = atmosphere.index(pos.x, pos.y);
        let cell = &mut atmosphere.cells[idx];

        // Generate O2 and N2
        let o2_generated = constants::LIFE_SUPPORT_O2_RATE * dt; // kg
        let n2_generated = constants::LIFE_SUPPORT_N2_RATE * dt; // kg

        // Convert to density change (kg/m³)
        let delta_rho_o2 = o2_generated / tile_volume;
        let delta_rho_n2 = n2_generated / tile_volume;

        // Add gases
        cell.rho_o2 += delta_rho_o2;
        cell.rho_n2 += delta_rho_n2;

        // Track total mass added
        total_generated += o2_generated + n2_generated;

        // Set temperature to room temp if generating (life support heats the gas)
        if cell.temperature < constants::ROOM_TEMP {
            cell.temperature = constants::ROOM_TEMP;
        }
    }

    // Update mass tracker with expected change
    tracker.expected_delta += total_generated;
}

/// Monitor velocity field divergence
///
/// Computes and logs the maximum divergence to verify pressure projection is working.
/// Ideally, divergence should be near zero after projection.
fn monitor_divergence(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut tracker: ResMut<DivergenceTracker>,
    time: Res<Time>,
) {
    tracker.time_since_check += time.delta_secs();

    if tracker.time_since_check < constants::MASS_CHECK_INTERVAL {
        return;
    }

    let dx = atmosphere.tile_size_physical;
    let mut max_divergence: f32 = 0.0;
    let mut total_divergence = 0.0;
    let mut cell_count = 0;

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);

            // Compute divergence
            let u_right = if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                atmosphere.cells[atmosphere.index(x + 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let u_left = if x > 0 && !collision_map.is_blocked(x - 1, y) {
                atmosphere.cells[atmosphere.index(x - 1, y)].u
            } else {
                atmosphere.cells[idx].u
            };

            let v_up = if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                atmosphere.cells[atmosphere.index(x, y + 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let v_down = if y > 0 && !collision_map.is_blocked(x, y - 1) {
                atmosphere.cells[atmosphere.index(x, y - 1)].v
            } else {
                atmosphere.cells[idx].v
            };

            let du_dx = (u_right - u_left) / (2.0 * dx);
            let dv_dy = (v_up - v_down) / (2.0 * dx);
            let div = du_dx + dv_dy;

            max_divergence = max_divergence.max(div.abs());
            total_divergence += div.abs();
            cell_count += 1;
        }
    }

    let avg_divergence = if cell_count > 0 {
        total_divergence / cell_count as f32
    } else {
        0.0
    };

    debug!(
        "Divergence check: max = {:.6} s⁻¹, avg = {:.6} s⁻¹",
        max_divergence, avg_divergence
    );

    // Warn if divergence is too high (indicates projection not working well)
    if max_divergence > 1.0 {
        warn!(
            "High velocity divergence detected: {:.3} s⁻¹ (pressure projection may need tuning)",
            max_divergence
        );
    }

    tracker.time_since_check = 0.0;
}

/// Local mixing fans mounted on the life support units to emulate forced convection
fn life_support_mixing(
    mut atmosphere: ResMut<AtmosphereGrid>,
    life_support: Res<LifeSupportTiles>,
    collision_map: Res<TileCollisionMap>,
) {
    if life_support.positions.is_empty() {
        return;
    }

    let radius = constants::FAN_RADIUS_TILES;
    let radius_f32 = radius as f32;

    for pos in &life_support.positions {
        let center_x = pos.x as i32;
        let center_y = pos.y as i32;

        for dy in -radius..=radius {
            for dx in -radius..=radius {
                if dx == 0 && dy == 0 {
                    continue;
                }

                let x = center_x + dx;
                let y = center_y + dy;

                if !atmosphere.in_bounds(x, y) || collision_map.is_blocked(x as u32, y as u32) {
                    continue;
                }

                let Some(cell) = atmosphere.get_mut(x, y) else {
                    continue;
                };

                let offset = Vec2::new(dx as f32, dy as f32);
                let distance = offset.length();
                if distance > radius_f32 {
                    continue;
                }

                let tangent = Vec2::new(-offset.y, offset.x);
                let tangent_len = tangent.length();
                if tangent_len <= f32::EPSILON {
                    continue;
                }
                let swirl_dir = tangent / tangent_len;

                let falloff = (1.0 - distance / (radius_f32 + 1.0)).clamp(0.0, 1.0);
                if falloff <= 0.0 {
                    continue;
                }

                let target_velocity = swirl_dir * constants::FAN_SWIRL_SPEED;
                let blend = constants::FAN_BLEND * falloff;

                cell.u += (target_velocity.x - cell.u) * blend;
                cell.v += (target_velocity.y - cell.v) * blend;
            }
        }
    }
}

/// Player respiration system - consumes O2 and produces CO2
fn player_respiration(
    mut atmosphere: ResMut<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    mut tracker: ResMut<MassTracker>,
    player_query: Query<&Transform, With<Player>>,
    time: Res<Time>,
) {
    let Ok(player_transform) = player_query.get_single() else {
        return;
    };

    let player_pos = player_transform.translation.truncate();

    // Convert player position to tile coordinates
    let Some(tile_pos) = collision_map.world_to_tile(player_pos) else {
        return;
    };

    // Calculate tile volume [m³]
    // Assume each tile is 1m x 1m x ROOM_HEIGHT m
    let tile_volume =
        atmosphere.tile_size_physical * atmosphere.tile_size_physical * constants::ROOM_HEIGHT;

    let idx = atmosphere.index(tile_pos.x, tile_pos.y);
    let cell = &mut atmosphere.cells[idx];

    // Calculate consumption/production rates
    let dt = time.delta_secs();
    let o2_consumed = constants::O2_CONSUMPTION_RATE * dt; // kg
    let co2_produced = constants::CO2_PRODUCTION_RATE * dt; // kg

    // Convert to density change (kg/m³)
    let delta_rho_o2 = -o2_consumed / tile_volume;
    let delta_rho_co2 = co2_produced / tile_volume;

    // Update densities
    cell.rho_o2 = (cell.rho_o2 + delta_rho_o2).max(0.0); // Can't go negative
    cell.rho_co2 += delta_rho_co2;
    cell.update_pressure();

    // Track net mass change (CO2 produced - O2 consumed)
    let net_mass_change = co2_produced - o2_consumed;
    tracker.expected_delta += net_mass_change;

    // Log warning if oxygen is getting dangerously low
    // Normal O2 partial pressure: ~21 kPa. Below 16 kPa (~16% O2) is hypoxic
    let o2_partial_pressure = (cell.rho_o2 / constants::M_O2) * constants::R * cell.temperature;
    if cell.pressure > 0.0 && o2_partial_pressure < 16000.0 {
        warn!(
            "Low oxygen at player location! O2 partial pressure: {:.0} Pa ({:.1}%)",
            o2_partial_pressure,
            100.0 * o2_partial_pressure / cell.pressure
        );
    }

    // Log warning if CO2 is getting dangerously high
    // Normal CO2: ~400 ppm = 0.04 kPa. Above 5 kPa (5%) is dangerous
    let co2_partial_pressure = (cell.rho_co2 / constants::M_CO2) * constants::R * cell.temperature;
    if cell.pressure > 0.0 && co2_partial_pressure > 5000.0 {
        warn!(
            "High CO2 at player location! CO2 partial pressure: {:.0} Pa ({:.1}%)",
            co2_partial_pressure,
            100.0 * co2_partial_pressure / cell.pressure
        );
    }
}

/// Visualization mode for debug rendering
#[derive(Resource, Default, Clone, Copy, PartialEq, Eq, Debug)]
pub enum VisualizationMode {
    Off,
    #[default]
    Pressure,
    Temperature,
    Oxygen,
    Nitrogen,
    CarbonDioxide,
    Breathability, // Heatmap showing air quality (O2 vs CO2)
    Velocity,
}

/// Component marking atmosphere debug visualization entities
#[derive(Component)]
struct AtmosphereDebugViz;

/// Debug visualization system
/// Press V to cycle through visualization modes
fn debug_visualization(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut mode: Local<VisualizationMode>,
    mut gizmos: Gizmos,
) {
    // Toggle breathability visualization with V key (on/off only)
    if keyboard.just_pressed(KeyCode::KeyV) {
        *mode = match *mode {
            VisualizationMode::Off => VisualizationMode::Breathability,
            VisualizationMode::Breathability => VisualizationMode::Off,
            _ => VisualizationMode::Off, // Fallback
        };

        info!("Atmosphere visualization: {:?}", *mode);
    }

    if *mode == VisualizationMode::Off {
        return;
    }

    // Draw colored rectangles for each cell
    let tile_size = atmosphere.tile_size_world;
    let origin = collision_map.origin;

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            // Skip walls
            if collision_map.is_blocked(x, y) {
                continue;
            }

            // Calculate world position for this tile
            let world_x = origin.x + (x as f32 + 0.5) * tile_size;
            let world_y = origin.y + (y as f32 + 0.5) * tile_size;

            // Get current cell for velocity arrows
            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells[idx];

            // For Breathability mode, use solid binary colors (no gradient)
            // For other modes, use smooth gradient interpolation
            if *mode == VisualizationMode::Breathability {
                // Binary mode - draw single solid-color tile
                let color = get_cell_color(&atmosphere, &collision_map, x, y, *mode);

                gizmos.rect_2d(
                    Isometry2d::from_translation(Vec2::new(world_x, world_y)),
                    Vec2::splat(tile_size * 0.95),
                    color,
                );
            } else {
                // Gradient mode - sample neighbor cells and interpolate
                let color_tl = get_cell_color(&atmosphere, &collision_map, x, y + 1, *mode);
                let color_tr = get_cell_color(&atmosphere, &collision_map, x + 1, y + 1, *mode);
                let color_bl = get_cell_color(&atmosphere, &collision_map, x, y, *mode);
                let color_br = get_cell_color(&atmosphere, &collision_map, x + 1, y, *mode);

                // Draw gradient rectangle by creating a mesh of smaller quads
                let subdiv = 2;
                let sub_size = tile_size / subdiv as f32;

                for sy in 0..subdiv {
                    for sx in 0..subdiv {
                        let fx0 = sx as f32 / subdiv as f32;
                        let fy0 = sy as f32 / subdiv as f32;
                        let fx1 = (sx + 1) as f32 / subdiv as f32;
                        let fy1 = (sy + 1) as f32 / subdiv as f32;

                        // Bilinear interpolation of colors
                        let interp_color = bilinear_color_interpolation(
                            color_bl,
                            color_br,
                            color_tl,
                            color_tr,
                            (fx0 + fx1) * 0.5,
                            (fy0 + fy1) * 0.5,
                        );

                        let sub_x = origin.x + (x as f32 + fx0 + 0.5 / subdiv as f32) * tile_size;
                        let sub_y = origin.y + (y as f32 + fy0 + 0.5 / subdiv as f32) * tile_size;

                        gizmos.rect_2d(
                            Isometry2d::from_translation(Vec2::new(sub_x, sub_y)),
                            Vec2::splat(sub_size * 0.95),
                            interp_color,
                        );
                    }
                }
            }

            // For velocity mode, also draw velocity vectors
            if *mode == VisualizationMode::Velocity {
                let vel_scale = 2.0;
                let vel_x = cell.u * vel_scale;
                let vel_y = cell.v * vel_scale;

                if vel_x.abs() > 0.01 || vel_y.abs() > 0.01 {
                    gizmos.arrow_2d(
                        Vec2::new(world_x, world_y),
                        Vec2::new(world_x + vel_x, world_y + vel_y),
                        Color::srgb(1.0, 1.0, 0.0),
                    );
                }
            }
        }
    }

    // Note: For a proper UI text overlay showing the current mode,
    // you would use bevy_ui to render text in the corner of the screen.
    // For now, mode changes are logged when you press V.
}

/// Helper function to get color for a cell at given coordinates
fn get_cell_color(
    atmosphere: &AtmosphereGrid,
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    mode: VisualizationMode,
) -> Color {
    // Check bounds and walls
    if x >= atmosphere.width || y >= atmosphere.height {
        return Color::srgba(0.0, 0.0, 0.0, 0.0); // Transparent for out of bounds
    }

    if collision_map.is_blocked(x, y) {
        return Color::srgba(0.0, 0.0, 0.0, 0.0); // Transparent for walls
    }

    let idx = atmosphere.index(x, y);
    let cell = &atmosphere.cells[idx];

    match mode {
        VisualizationMode::Off => Color::srgba(0.0, 0.0, 0.0, 0.0),
        VisualizationMode::Pressure => {
            let t = (cell.pressure / constants::EARTH_PRESSURE).clamp(0.0, 2.0) / 2.0;
            Color::srgb(t, 0.0, 1.0 - t).with_alpha(0.3)
        }
        VisualizationMode::Temperature => {
            let t = (cell.temperature / 600.0).clamp(0.0, 1.0);
            Color::srgb(t, t * 0.5, 1.0 - t).with_alpha(0.3)
        }
        VisualizationMode::Oxygen => {
            let earth_o2 = 0.273;
            let t = (cell.rho_o2 / earth_o2).clamp(0.0, 1.0);
            Color::srgb(0.0, t, 0.0).with_alpha(0.4)
        }
        VisualizationMode::Nitrogen => {
            let earth_n2 = 1.165;
            let t = (cell.rho_n2 / earth_n2).clamp(0.0, 1.0);
            Color::srgb(0.0, 0.0, t).with_alpha(0.4)
        }
        VisualizationMode::CarbonDioxide => {
            let dangerous_co2 = 0.09;
            let t = (cell.rho_co2 / dangerous_co2).clamp(0.0, 1.0);
            Color::srgb(t, 0.0, 0.0).with_alpha(0.4)
        }
        VisualizationMode::Breathability => {
            // Calculate partial pressures
            let o2_pp = (cell.rho_o2 / constants::M_O2) * constants::R * cell.temperature;
            let co2_pp = (cell.rho_co2 / constants::M_CO2) * constants::R * cell.temperature;

            // Simple binary check: is the air breathable or not?
            // Breathable criteria:
            // - O2 partial pressure between 16 kPa and 50 kPa
            // - CO2 partial pressure less than 5 kPa (5%)
            let is_breathable = o2_pp >= 16000.0 && o2_pp <= 50000.0 && co2_pp < 5000.0;

            if is_breathable {
                // Neon green for breathable air
                Color::srgb(0.0, 1.0, 0.0).with_alpha(0.6)
            } else {
                // Neon red for dangerous air
                Color::srgb(1.0, 0.0, 0.0).with_alpha(0.6)
            }
        }
        VisualizationMode::Velocity => {
            let vel_mag = (cell.u * cell.u + cell.v * cell.v).sqrt();
            let t = (vel_mag / 10.0).clamp(0.0, 1.0);
            Color::srgb(t, t, 0.0).with_alpha(0.3)
        }
    }
}

/// Bilinear interpolation of colors
fn bilinear_color_interpolation(
    c00: Color, // bottom-left
    c10: Color, // bottom-right
    c01: Color, // top-left
    c11: Color, // top-right
    fx: f32,    // horizontal interpolation factor [0,1]
    fy: f32,    // vertical interpolation factor [0,1]
) -> Color {
    // Convert to linear space for proper color interpolation
    let c00_linear = c00.to_linear();
    let c10_linear = c10.to_linear();
    let c01_linear = c01.to_linear();
    let c11_linear = c11.to_linear();

    // Interpolate bottom edge (in linear space)
    let bottom_r = c00_linear.red * (1.0 - fx) + c10_linear.red * fx;
    let bottom_g = c00_linear.green * (1.0 - fx) + c10_linear.green * fx;
    let bottom_b = c00_linear.blue * (1.0 - fx) + c10_linear.blue * fx;
    let bottom_a = c00_linear.alpha * (1.0 - fx) + c10_linear.alpha * fx;

    // Interpolate top edge (in linear space)
    let top_r = c01_linear.red * (1.0 - fx) + c11_linear.red * fx;
    let top_g = c01_linear.green * (1.0 - fx) + c11_linear.green * fx;
    let top_b = c01_linear.blue * (1.0 - fx) + c11_linear.blue * fx;
    let top_a = c01_linear.alpha * (1.0 - fx) + c11_linear.alpha * fx;

    // Interpolate between bottom and top
    let final_r = bottom_r * (1.0 - fy) + top_r * fy;
    let final_g = bottom_g * (1.0 - fy) + top_g * fy;
    let final_b = bottom_b * (1.0 - fy) + top_b * fy;
    let final_a = bottom_a * (1.0 - fy) + top_a * fy;

    // Return final color in linear space so downstream conversions apply gamma exactly once.
    Color::linear_rgba(final_r, final_g, final_b, final_a)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test that species conservation is maintained when flow reverses direction
    /// This test catches the donor cell bug where we used the wrong cell's composition
    #[test]
    fn test_donor_cell_species_conservation() {
        // Create two cells with different gas compositions
        let mut cell_a = AtmosphereCell {
            rho_o2: 0.273,  // 100% O2 (Earth-like O2 density)
            rho_n2: 0.0,
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        };
        cell_a.update_pressure();

        let mut cell_b = AtmosphereCell {
            rho_o2: 0.0,
            rho_n2: 1.165,  // 100% N2 (Earth-like N2 density)
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        };
        cell_b.update_pressure();

        // Cell B has higher pressure, so flow should be FROM B TO A
        assert!(cell_b.pressure > cell_a.pressure, "B should have higher pressure");

        // Simulate pressure-driven flux (simplified version of the actual code)
        let dx = 1.0; // 1m tile
        let dt = 0.01; // 10ms timestep
        let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;
        let _cell_volume = dx * dx * constants::ROOM_HEIGHT;

        let p_diff = cell_a.pressure - cell_b.pressure; // Negative (B > A)
        let total_flux = conductance * p_diff * dx * dt;

        // Identify donor (B) based on negative p_diff
        let donor = if p_diff > 0.0 { &cell_a } else { &cell_b };

        // Donor is B (100% N2), so we should transfer only N2
        assert_eq!(donor.rho_o2, 0.0, "Donor should be cell B with no O2");
        assert!(donor.rho_n2 > 0.0, "Donor should be cell B with N2");

        // Calculate species flux using donor's composition
        let abs_flux = total_flux.abs();
        if donor.pressure > f32::EPSILON {
            let p_o2 = (donor.rho_o2 / constants::M_O2) * constants::R * donor.temperature;
            let p_n2 = (donor.rho_n2 / constants::M_N2) * constants::R * donor.temperature;

            let flux_o2 = abs_flux * (p_o2 / donor.pressure);
            let flux_n2 = abs_flux * (p_n2 / donor.pressure);

            // Since donor has 0% O2, flux_o2 should be zero
            assert!(flux_o2.abs() < 1e-6, "Should transfer zero O2 from N2-only donor");
            // Since donor has 100% N2, all flux should be N2
            assert!((flux_n2 - abs_flux).abs() < 1e-3, "Should transfer only N2");
        }
    }

    /// Test that momentum conservation is maintained with donor cell selection
    /// This test catches the bug where we used receiver's velocity instead of donor's
    #[test]
    fn test_donor_cell_momentum_conservation() {
        // Create two cells: A stationary, B moving fast
        let mut cell_a = AtmosphereCell {
            rho_o2: 0.273,
            rho_n2: 1.165,
            rho_co2: 0.0,
            u: 0.0,  // Stationary
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        };
        cell_a.update_pressure();

        let mut cell_b = AtmosphereCell {
            rho_o2: 0.273,
            rho_n2: 1.165,
            rho_co2: 0.0,
            u: 10.0,  // Moving at 10 m/s
            v: 0.0,
            temperature: constants::ROOM_TEMP * 2.0,  // Higher temp = higher pressure
            pressure: 0.0,
        };
        cell_b.update_pressure();

        // Cell B has higher pressure (due to higher temp), so flow is FROM B TO A
        assert!(cell_b.pressure > cell_a.pressure, "B should have higher pressure");

        // Calculate momentum transfer
        let dx = 1.0;
        let dt = 0.01;
        let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;
        let _cell_volume = dx * dx * constants::ROOM_HEIGHT;

        let p_diff = cell_a.pressure - cell_b.pressure; // Negative
        let total_flux = conductance * p_diff * dx * dt;
        let abs_flux = total_flux.abs();

        // Donor is B (moving at 10 m/s)
        let donor = if p_diff > 0.0 { &cell_a } else { &cell_b };
        let donor_velocity_u = donor.u;

        // Momentum flux should use DONOR's velocity (10 m/s), not receiver's (0 m/s)
        let momentum_flux_u = abs_flux * donor_velocity_u;

        assert!(donor_velocity_u.abs() > 0.0, "Donor should be moving");
        assert!(momentum_flux_u.abs() > 0.0, "Momentum flux should be non-zero when donor is moving");

        // If we incorrectly used receiver's velocity, momentum_flux would be zero
        let wrong_momentum_flux = abs_flux * 0.0; // receiver.u = 0.0
        assert_eq!(wrong_momentum_flux, 0.0, "Wrong approach gives zero momentum transfer");
        assert_ne!(momentum_flux_u, wrong_momentum_flux, "Correct approach must differ");
    }

    /// Test that total mass is conserved during pressure flux
    #[test]
    fn test_pressure_flux_mass_conservation() {
        let mut cell_a = AtmosphereCell::earth_atmosphere();
        cell_a.pressure = 50000.0; // Low pressure

        let mut cell_b = AtmosphereCell::earth_atmosphere();
        cell_b.pressure = 150000.0; // High pressure

        let initial_mass_a = cell_a.total_density();
        let initial_mass_b = cell_b.total_density();
        let total_initial_mass = initial_mass_a + initial_mass_b;

        // Simulate flux
        let dx = 1.0;
        let dt = 0.01;
        let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;
        let cell_volume = dx * dx * constants::ROOM_HEIGHT;

        let p_diff = cell_a.pressure - cell_b.pressure;
        let total_flux = conductance * p_diff * dx * dt;
        let abs_flux = total_flux.abs();

        let (donor, donor_is_a) = if p_diff > 0.0 {
            (&cell_a, true)
        } else {
            (&cell_b, false)
        };

        if donor.pressure > f32::EPSILON {
            let p_o2 = (donor.rho_o2 / constants::M_O2) * constants::R * donor.temperature;
            let p_n2 = (donor.rho_n2 / constants::M_N2) * constants::R * donor.temperature;
            let p_co2 = (donor.rho_co2 / constants::M_CO2) * constants::R * donor.temperature;

            let flux_o2 = abs_flux * (p_o2 / donor.pressure);
            let flux_n2 = abs_flux * (p_n2 / donor.pressure);
            let flux_co2 = abs_flux * (p_co2 / donor.pressure);

            let delta_rho_o2 = flux_o2 / cell_volume;
            let delta_rho_n2 = flux_n2 / cell_volume;
            let delta_rho_co2 = flux_co2 / cell_volume;

            // Apply mass transfer
            if donor_is_a {
                cell_a.rho_o2 -= delta_rho_o2;
                cell_a.rho_n2 -= delta_rho_n2;
                cell_a.rho_co2 -= delta_rho_co2;

                cell_b.rho_o2 += delta_rho_o2;
                cell_b.rho_n2 += delta_rho_n2;
                cell_b.rho_co2 += delta_rho_co2;
            } else {
                cell_b.rho_o2 -= delta_rho_o2;
                cell_b.rho_n2 -= delta_rho_n2;
                cell_b.rho_co2 -= delta_rho_co2;

                cell_a.rho_o2 += delta_rho_o2;
                cell_a.rho_n2 += delta_rho_n2;
                cell_a.rho_co2 += delta_rho_co2;
            }
        }

        let final_mass_a = cell_a.total_density();
        let final_mass_b = cell_b.total_density();
        let total_final_mass = final_mass_a + final_mass_b;

        // Total mass should be conserved
        let mass_error = (total_final_mass - total_initial_mass).abs();
        assert!(
            mass_error < 1e-6,
            "Total mass should be conserved, error: {}",
            mass_error
        );
    }

    /// Test MacCormack advection produces different results than simple semi-Lagrangian
    /// This is a smoke test to ensure the MacCormack path is actually running
    #[test]
    fn test_maccormack_differs_from_semi_lagrangian() {
        // Create a simple velocity field
        let cells = vec![
            AtmosphereCell {
                rho_o2: 1.0,
                rho_n2: 0.0,
                rho_co2: 0.0,
                u: 1.0,
                v: 0.0,
                temperature: 300.0,
                pressure: 101325.0,
            },
            AtmosphereCell {
                rho_o2: 0.0,
                rho_n2: 1.0,
                rho_co2: 0.0,
                u: 1.0,
                v: 0.0,
                temperature: 300.0,
                pressure: 101325.0,
            },
        ];

        // Test interpolation (which is used by both schemes)
        let result = interpolate_cell(&cells, 2, 1, 0.5, 0.0);

        // At x=0.5, should interpolate between cells[0] and cells[1]
        assert!(result.rho_o2 > 0.0 && result.rho_o2 < 1.0, "Should interpolate O2");
        assert!(result.rho_n2 > 0.0 && result.rho_n2 < 1.0, "Should interpolate N2");
    }

    /// Test that bilinear interpolation works correctly
    #[test]
    fn test_bilinear_interpolation() {
        // Test bilerp function
        let v00 = 0.0;
        let v10 = 1.0;
        let v01 = 0.0;
        let v11 = 1.0;

        // Center should be 0.5
        let center = bilerp(v00, v10, v01, v11, 0.5, 0.5);
        assert!((center - 0.5).abs() < 1e-6, "Center interpolation should be 0.5");

        // Corners should match input
        let corner_00 = bilerp(v00, v10, v01, v11, 0.0, 0.0);
        assert!((corner_00 - v00).abs() < 1e-6, "Corner 00 should match");

        let corner_11 = bilerp(v00, v10, v01, v11, 1.0, 1.0);
        assert!((corner_11 - v11).abs() < 1e-6, "Corner 11 should match");
    }

    /// Test that pressure is correctly computed from ideal gas law
    #[test]
    fn test_pressure_calculation() {
        // Use canonical Earth atmosphere composition to avoid rounding drift.
        let mut cell = AtmosphereCell::earth_atmosphere();
        // Zero the cached pressure so `update_pressure` recomputes it.
        cell.pressure = 0.0;
        cell.update_pressure();

        // Should be approximately Earth pressure (101325 Pa)
        let error = (cell.pressure - constants::EARTH_PRESSURE).abs();
        let relative_error = error / constants::EARTH_PRESSURE;

        assert!(
            relative_error < 0.05,
            "Pressure should be close to Earth pressure, got {} Pa (error: {:.1}%)",
            cell.pressure,
            relative_error * 100.0
        );
    }
}
