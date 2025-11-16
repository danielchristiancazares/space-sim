use crate::atmosphere::constants;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::Resource;

/// Atmospheric state for a single grid cell.
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
    #[allow(dead_code)]
    /// Create cell with Earth-like atmosphere.
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

    /// Create vacuum cell.
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

    /// Compute pressure from current state using ideal gas law.
    /// p = (ρ_O2/M_O2 + ρ_N2/M_N2 + ρ_CO2/M_CO2) * R * T
    #[inline]
    pub fn compute_pressure(&self) -> f32 {
        let n_o2 = self.rho_o2 / constants::M_O2;
        let n_n2 = self.rho_n2 / constants::M_N2;
        let n_co2 = self.rho_co2 / constants::M_CO2;
        (n_o2 + n_n2 + n_co2) * constants::R * self.temperature
    }

    /// Update cached pressure value.
    #[inline]
    pub fn update_pressure(&mut self) {
        self.pressure = self.compute_pressure();
    }
}

/// Grid-based atmospheric simulation resource.
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

    #[allow(dead_code)]
    /// Initialize with Earth-like atmosphere in open areas.
    pub fn initialize_earth_atmosphere(&mut self, collision_map: &TileCollisionMap) {
        for y in 0..self.height {
            for x in 0..self.width {
                let idx = self.index(x, y);

                if !collision_map.is_blocked(x, y) {
                    self.cells[idx] = AtmosphereCell::earth_atmosphere();
                } else {
                    self.cells[idx] = AtmosphereCell::vacuum();
                }
            }
        }
    }

    /// Initialize entire room as vacuum.
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

    #[allow(dead_code)]
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pressure_calculation() {
        let mut cell = AtmosphereCell::earth_atmosphere();
        cell.pressure = 0.0;
        cell.update_pressure();

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
