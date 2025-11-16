#![allow(dead_code)]

/// Physical constants and tuning parameters for the atmospheric simulation.
///
/// This module is the single source of truth for all simulation constants so that
/// other modules can reference them without duplicating values.
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

/// Effective bulk gas diffusivity used for scalar mixing [m²/s]
pub const GAS_DIFFUSIVITY: f32 = 5.0e-4;

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

/// Earth atmosphere O2 density at sea level [kg/m³]
pub const EARTH_RHO_O2: f32 = (EARTH_PRESSURE * EARTH_O2_FRACTION * M_O2) / (R * ROOM_TEMP);

/// Earth atmosphere N2 density at sea level [kg/m³]
pub const EARTH_RHO_N2: f32 = (EARTH_PRESSURE * EARTH_N2_FRACTION * M_N2) / (R * ROOM_TEMP);

/// Earth atmosphere CO2 density at sea level [kg/m³]
pub const EARTH_RHO_CO2: f32 = (EARTH_PRESSURE * EARTH_CO2_FRACTION * M_CO2) / (R * ROOM_TEMP);

/// CO2 density corresponding to a 5 kPa partial pressure threshold [kg/m³]
pub const CO2_DANGER_DENSITY: f32 = (5000.0 * M_CO2) / (R * ROOM_TEMP);

/// Kinematic viscosity of air at ~20°C [m²/s] (ν)
pub const KINEMATIC_VISCOSITY: f32 = 1.5e-5;

/// Thermal diffusivity of air at ~20°C [m²/s] (α_th)
pub const THERMAL_DIFFUSIVITY: f32 = 2.1e-5;

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
/// 0.35 was chosen so that air parcels within the fan radius reach ~90% of the
/// target swirl in roughly three simulation frames (1 - 0.35^3 ≈ 0.96), which
/// matched playtest expectations for “noticeable but not instant” circulation.
pub const FAN_BLEND: f32 = 0.35;

/// Mass conservation check interval (seconds)
/// How often to verify total mass conservation
pub const MASS_CHECK_INTERVAL: f32 = 5.0;

/// Divergence monitoring interval (seconds)
pub const DIVERGENCE_CHECK_INTERVAL: f32 = 5.0;

/// Mass conservation tolerance (relative error)
/// Warn if mass changes by more than this fraction (excluding life support/respiration)
pub const MASS_TOLERANCE: f32 = 0.01; // 1%

/// CFL number for timestep stability
/// CFL = (u*dt)/dx should be < 1 for explicit schemes
/// Using 0.5 provides good stability margin
pub const CFL_NUMBER: f32 = 0.5;

/// Ratio of specific heats (Cp/Cv) for air [dimensionless]
/// For diatomic gases (O₂, N₂): γ ≈ 1.4
/// Used in compression/expansion heating: ∂T/∂t = -(γ-1)T∇·u
pub const GAMMA: f32 = 1.4;

/// Specific gas constant for air [J/(kg·K)]
/// R_specific = R_universal / M_air where M_air ≈ 0.029 kg/mol
/// Used in sound speed calculation: c = √(γRT)
pub const R_SPECIFIC_AIR: f32 = 287.05;

/// CFL diagnostic mode (safety/paranoid mode for debugging)
/// Use this if standard CFL_NUMBER = 0.5 shows instability
pub const CFL_DIAGNOSTIC: f32 = 0.3;

/// Cosmic microwave background temperature [K]
/// Absolute physical minimum for temperature clamping
pub const T_CMB: f32 = 2.7;

// ============================================================================
// Artificial Viscosity Presets for Shock Regularization
// ============================================================================
// At coarse grid resolution (1m cells), molecular viscosity is insufficient
// to regularize supersonic expansion shocks. Artificial viscosity represents
// unresolved subgrid turbulence and prevents numerical oscillations.
//
// Formula: ν_artificial = COEFF × Δx × |u_max|
//
// Quality modes:
//   - QUALITY (0.05):     Minimal damping, preserves fine structures
//   - PRODUCTION (0.1):   Balanced damping (DEFAULT for research)
//   - SAFETY (0.2):       Strong damping, very stable

/// Artificial viscosity coefficient: minimal damping
pub const ARTIFICIAL_VISCOSITY_QUALITY: f32 = 0.05;

/// Artificial viscosity coefficient: balanced damping (DEFAULT)
pub const ARTIFICIAL_VISCOSITY_PRODUCTION: f32 = 0.1;

/// Artificial viscosity coefficient: strong damping
pub const ARTIFICIAL_VISCOSITY_SAFETY: f32 = 0.2;

/// Active artificial viscosity setting
pub const ARTIFICIAL_VISCOSITY: f32 = ARTIFICIAL_VISCOSITY_PRODUCTION;
