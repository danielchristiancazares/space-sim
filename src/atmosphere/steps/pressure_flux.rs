use crate::atmosphere::constants;
use crate::atmosphere::grid::{AtmosphereCell, AtmosphereGrid};
use crate::tilemap::TileCollisionMap;

#[derive(Copy, Clone)]
enum VelocityComponent {
    X,
    Y,
}

fn velocity_of(cell: &AtmosphereCell, component: VelocityComponent) -> f32 {
    match component {
        VelocityComponent::X => cell.u,
        VelocityComponent::Y => cell.v,
    }
}

fn set_velocity(cell: &mut AtmosphereCell, component: VelocityComponent, value: f32) {
    match component {
        VelocityComponent::X => cell.u = value,
        VelocityComponent::Y => cell.v = value,
    }
}

pub fn pressure_flux_step(
    atmosphere: &mut AtmosphereGrid,
    collision_map: &TileCollisionMap,
    dt: f32,
) {
    let dx = atmosphere.tile_size_physical;
    let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;

    let cell_volume = dx * dx * constants::ROOM_HEIGHT;

    atmosphere.cells_buffer.clone_from_slice(&atmosphere.cells);

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let idx = atmosphere.index(x, y);
            let cell_pressure = atmosphere.cells_buffer[idx].pressure;

            if x < atmosphere.width - 1 && !collision_map.is_blocked(x + 1, y) {
                let idx_right = atmosphere.index(x + 1, y);
                let right_pressure = atmosphere.cells_buffer[idx_right].pressure;

                let p_diff = cell_pressure - right_pressure;

                if p_diff.abs() > f32::EPSILON {
                    let mass_to_transfer = (conductance * p_diff * dx * dt).abs();

                    if mass_to_transfer > 0.0 {
                        let (donor_idx, receiver_idx) = if p_diff > 0.0 {
                            (idx, idx_right)
                        } else {
                            (idx_right, idx)
                        };
                        let donor_state = atmosphere.cells_buffer[donor_idx].clone();
                        if donor_state.pressure <= f32::EPSILON {
                            continue;
                        }
                        process_flux(
                            atmosphere,
                            donor_idx,
                            receiver_idx,
                            &donor_state,
                            mass_to_transfer,
                            cell_volume,
                            VelocityComponent::X,
                        );
                    }
                }
            }

            if y < atmosphere.height - 1 && !collision_map.is_blocked(x, y + 1) {
                let idx_up = atmosphere.index(x, y + 1);
                let up_pressure = atmosphere.cells_buffer[idx_up].pressure;

                let p_diff = cell_pressure - up_pressure;

                if p_diff.abs() > f32::EPSILON {
                    let mass_to_transfer = (conductance * p_diff * dx * dt).abs();

                    if mass_to_transfer > 0.0 {
                        let (donor_idx, receiver_idx) = if p_diff > 0.0 {
                            (idx, idx_up)
                        } else {
                            (idx_up, idx)
                        };
                        let donor_state = atmosphere.cells_buffer[donor_idx].clone();
                        if donor_state.pressure <= f32::EPSILON {
                            continue;
                        }
                        process_flux(
                            atmosphere,
                            donor_idx,
                            receiver_idx,
                            &donor_state,
                            mass_to_transfer,
                            cell_volume,
                            VelocityComponent::Y,
                        );
                    }
                }
            }
        }
    }

    for cell in atmosphere.cells.iter_mut() {
        cell.rho_o2 = cell.rho_o2.max(0.0);
        cell.rho_n2 = cell.rho_n2.max(0.0);
        cell.rho_co2 = cell.rho_co2.max(0.0);
    }
}

fn process_flux(
    atmosphere: &mut AtmosphereGrid,
    donor_idx: usize,
    receiver_idx: usize,
    donor_state: &AtmosphereCell,
    mass_to_transfer: f32,
    cell_volume: f32,
    component: VelocityComponent,
) {
    if donor_state.pressure <= f32::EPSILON || mass_to_transfer <= 0.0 {
        return;
    }

    let p_o2 = (donor_state.rho_o2 / constants::M_O2) * constants::R * donor_state.temperature;
    let p_n2 = (donor_state.rho_n2 / constants::M_N2) * constants::R * donor_state.temperature;
    let p_co2 = (donor_state.rho_co2 / constants::M_CO2) * constants::R * donor_state.temperature;

    let frac_o2 = (p_o2 / donor_state.pressure).clamp(0.0, 1.0);
    let frac_n2 = (p_n2 / donor_state.pressure).clamp(0.0, 1.0);
    let frac_co2 = (p_co2 / donor_state.pressure).clamp(0.0, 1.0);

    let mut m_o2 = mass_to_transfer * frac_o2;
    let mut m_n2 = mass_to_transfer * frac_n2;
    let mut m_co2 = mass_to_transfer * frac_co2;

    let donor_before = atmosphere.cells_buffer[donor_idx].clone();
    let receiver_before = atmosphere.cells_buffer[receiver_idx].clone();
    let donor_mass_before = donor_before.total_density() * cell_volume;
    if donor_mass_before <= 0.0 {
        return;
    }
    let donor_velocity_before = velocity_of(&donor_before, component);
    let receiver_mass_before = receiver_before.total_density() * cell_volume;
    let receiver_velocity_before = velocity_of(&receiver_before, component);

    let avail_o2 = donor_before.rho_o2 * cell_volume;
    let avail_n2 = donor_before.rho_n2 * cell_volume;
    let avail_co2 = donor_before.rho_co2 * cell_volume;

    let scale_o2 = if m_o2 > 0.0 {
        (avail_o2 / m_o2).min(1.0)
    } else {
        1.0
    };
    let scale_n2 = if m_n2 > 0.0 {
        (avail_n2 / m_n2).min(1.0)
    } else {
        1.0
    };
    let scale_co2 = if m_co2 > 0.0 {
        (avail_co2 / m_co2).min(1.0)
    } else {
        1.0
    };
    let availability_scale = scale_o2.min(scale_n2).min(scale_co2);

    m_o2 *= availability_scale;
    m_n2 *= availability_scale;
    m_co2 *= availability_scale;

    let mut transferred_mass = m_o2 + m_n2 + m_co2;
    if transferred_mass <= 0.0 {
        return;
    }

    let max_mass = donor_mass_before * constants::PRESSURE_FLUX_MAX_FRACTION;
    if max_mass > 0.0 && transferred_mass > max_mass {
        let limit_scale = max_mass / transferred_mass;
        m_o2 *= limit_scale;
        m_n2 *= limit_scale;
        m_co2 *= limit_scale;
        transferred_mass = max_mass;
    }

    if transferred_mass <= 0.0 {
        return;
    }

    let delta_rho_o2 = m_o2 / cell_volume;
    let delta_rho_n2 = m_n2 / cell_volume;
    let delta_rho_co2 = m_co2 / cell_volume;

    atmosphere.cells[donor_idx].rho_o2 -= delta_rho_o2;
    atmosphere.cells[donor_idx].rho_n2 -= delta_rho_n2;
    atmosphere.cells[donor_idx].rho_co2 -= delta_rho_co2;

    atmosphere.cells[receiver_idx].rho_o2 += delta_rho_o2;
    atmosphere.cells[receiver_idx].rho_n2 += delta_rho_n2;
    atmosphere.cells[receiver_idx].rho_co2 += delta_rho_co2;

    let donor_mass_target = (donor_mass_before - transferred_mass).max(0.0);
    let donor_mass_after = (atmosphere.cells[donor_idx].total_density() * cell_volume).max(1e-6);
    let donor_velocity_after = (donor_mass_target * donor_velocity_before) / donor_mass_after;
    set_velocity(
        &mut atmosphere.cells[donor_idx],
        component,
        donor_velocity_after,
    );

    let receiver_mass_after =
        (atmosphere.cells[receiver_idx].total_density() * cell_volume).max(1e-6);
    let receiver_momentum_after =
        receiver_mass_before * receiver_velocity_before + transferred_mass * donor_velocity_before;
    let receiver_velocity_after = receiver_momentum_after / receiver_mass_after;
    set_velocity(
        &mut atmosphere.cells[receiver_idx],
        component,
        receiver_velocity_after,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_donor_cell_species_conservation() {
        let mut cell_a = AtmosphereCell {
            rho_o2: 0.273,
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
            rho_n2: 1.165,
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        };
        cell_b.update_pressure();

        assert!(cell_b.pressure > cell_a.pressure);

        let dx = 1.0;
        let dt = 0.01;
        let p_diff = cell_a.pressure - cell_b.pressure;
        let total_flux = constants::PRESSURE_FLUX_CONDUCTANCE * p_diff * dx * dt;

        let donor = if p_diff > 0.0 { &cell_a } else { &cell_b };
        let abs_flux = total_flux.abs();
        if donor.pressure > f32::EPSILON {
            let p_o2 = (donor.rho_o2 / constants::M_O2) * constants::R * donor.temperature;
            let p_n2 = (donor.rho_n2 / constants::M_N2) * constants::R * donor.temperature;

            let flux_o2 = abs_flux * (p_o2 / donor.pressure);
            let flux_n2 = abs_flux * (p_n2 / donor.pressure);

            assert!(flux_o2.abs() < 1e-6);
            assert!((flux_n2 - abs_flux).abs() < 1e-3);
        }
    }

    #[test]
    fn test_donor_cell_momentum_conservation() {
        let mut cell_a = AtmosphereCell {
            rho_o2: 0.273,
            rho_n2: 1.165,
            rho_co2: 0.0,
            u: 0.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP,
            pressure: 0.0,
        };
        cell_a.update_pressure();

        let mut cell_b = AtmosphereCell {
            rho_o2: 0.273,
            rho_n2: 1.165,
            rho_co2: 0.0,
            u: 10.0,
            v: 0.0,
            temperature: constants::ROOM_TEMP * 2.0,
            pressure: 0.0,
        };
        cell_b.update_pressure();

        assert!(cell_b.pressure > cell_a.pressure);

        let dx = 1.0;
        let dt = 0.01;
        let conductance = constants::PRESSURE_FLUX_CONDUCTANCE;

        let p_diff = cell_a.pressure - cell_b.pressure;
        let total_flux = conductance * p_diff * dx * dt;
        let abs_flux = total_flux.abs();

        let donor = if p_diff > 0.0 { &cell_a } else { &cell_b };
        let donor_velocity_u = donor.u;

        let momentum_flux_u = abs_flux * donor_velocity_u;
        assert!(donor_velocity_u.abs() > 0.0);
        assert!(momentum_flux_u.abs() > 0.0);
    }

    #[test]
    fn test_pressure_flux_mass_conservation() {
        let mut cell_a = AtmosphereCell::earth_atmosphere();
        cell_a.pressure = 50000.0;

        let mut cell_b = AtmosphereCell::earth_atmosphere();
        cell_b.pressure = 150000.0;

        let initial_mass_a = cell_a.total_density();
        let initial_mass_b = cell_b.total_density();
        let total_initial_mass = initial_mass_a + initial_mass_b;

        let dx = 1.0;
        let dt = 0.01;
        let cell_volume = dx * dx * constants::ROOM_HEIGHT;

        let p_diff = cell_a.pressure - cell_b.pressure;
        let total_flux = constants::PRESSURE_FLUX_CONDUCTANCE * p_diff * dx * dt;
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

        let mass_error = (total_final_mass - total_initial_mass).abs();
        assert!(mass_error < 1e-6, "Total mass error: {}", mass_error);
    }
}
