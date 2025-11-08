use crate::atmosphere::constants;
use crate::atmosphere::grid::AtmosphereGrid;
use crate::tilemap::TileCollisionMap;
use bevy::prelude::*;

#[allow(dead_code)]
#[derive(Resource, Default, Clone, Copy, PartialEq, Eq, Debug)]
pub enum VisualizationMode {
    Off,
    #[default]
    Pressure,
    Temperature,
    Oxygen,
    Nitrogen,
    CarbonDioxide,
    Breathability,
    Velocity,
}

pub fn debug_visualization(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut mode: Local<VisualizationMode>,
    mut gizmos: Gizmos,
) {
    if keyboard.just_pressed(KeyCode::KeyV) {
        *mode = match *mode {
            VisualizationMode::Off => VisualizationMode::Breathability,
            VisualizationMode::Breathability => VisualizationMode::Off,
            _ => VisualizationMode::Off,
        };

        info!("Atmosphere visualization: {:?}", *mode);
    }

    if *mode == VisualizationMode::Off {
        return;
    }

    let tile_size = atmosphere.tile_size_world;
    let origin = collision_map.origin;

    for y in 0..atmosphere.height {
        for x in 0..atmosphere.width {
            if collision_map.is_blocked(x, y) {
                continue;
            }

            let world_x = origin.x + (x as f32 + 0.5) * tile_size;
            let world_y = origin.y + (y as f32 + 0.5) * tile_size;
            let idx = atmosphere.index(x, y);
            let cell = &atmosphere.cells[idx];

            if *mode == VisualizationMode::Breathability {
                let color = get_cell_color(&atmosphere, &collision_map, x, y, *mode);
                gizmos.rect_2d(
                    Isometry2d::from_translation(Vec2::new(world_x, world_y)),
                    Vec2::splat(tile_size * 0.95),
                    color,
                );
            } else {
                let color_tl = get_cell_color(&atmosphere, &collision_map, x, y + 1, *mode);
                let color_tr = get_cell_color(&atmosphere, &collision_map, x + 1, y + 1, *mode);
                let color_bl = get_cell_color(&atmosphere, &collision_map, x, y, *mode);
                let color_br = get_cell_color(&atmosphere, &collision_map, x + 1, y, *mode);

                let subdiv = 2;
                let sub_size = tile_size / subdiv as f32;

                for sy in 0..subdiv {
                    for sx in 0..subdiv {
                        let fx0 = sx as f32 / subdiv as f32;
                        let fy0 = sy as f32 / subdiv as f32;
                        let fx1 = (sx + 1) as f32 / subdiv as f32;
                        let fy1 = (sy + 1) as f32 / subdiv as f32;

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
}

fn get_cell_color(
    atmosphere: &AtmosphereGrid,
    collision_map: &TileCollisionMap,
    x: u32,
    y: u32,
    mode: VisualizationMode,
) -> Color {
    if x >= atmosphere.width || y >= atmosphere.height {
        return Color::srgba(0.0, 0.0, 0.0, 0.0);
    }

    if collision_map.is_blocked(x, y) {
        return Color::srgba(0.0, 0.0, 0.0, 0.0);
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
            let t = (cell.rho_o2 / constants::EARTH_RHO_O2).clamp(0.0, 1.0);
            Color::srgb(0.0, t, 0.0).with_alpha(0.4)
        }
        VisualizationMode::Nitrogen => {
            let t = (cell.rho_n2 / constants::EARTH_RHO_N2).clamp(0.0, 1.0);
            Color::srgb(0.0, 0.0, t).with_alpha(0.4)
        }
        VisualizationMode::CarbonDioxide => {
            let t = (cell.rho_co2 / constants::CO2_DANGER_DENSITY).clamp(0.0, 1.0);
            Color::srgb(t, 0.0, 0.0).with_alpha(0.4)
        }
        VisualizationMode::Breathability => {
            let o2_pp = (cell.rho_o2 / constants::M_O2) * constants::R * cell.temperature;
            let co2_pp = (cell.rho_co2 / constants::M_CO2) * constants::R * cell.temperature;

            let is_breathable = o2_pp >= 16000.0 && o2_pp <= 50000.0 && co2_pp < 5000.0;

            if is_breathable {
                Color::srgb(0.0, 1.0, 0.0).with_alpha(0.6)
            } else {
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

fn bilinear_color_interpolation(
    c00: Color,
    c10: Color,
    c01: Color,
    c11: Color,
    fx: f32,
    fy: f32,
) -> Color {
    let c00_linear = c00.to_linear();
    let c10_linear = c10.to_linear();
    let c01_linear = c01.to_linear();
    let c11_linear = c11.to_linear();

    let bottom_r = c00_linear.red * (1.0 - fx) + c10_linear.red * fx;
    let bottom_g = c00_linear.green * (1.0 - fx) + c10_linear.green * fx;
    let bottom_b = c00_linear.blue * (1.0 - fx) + c10_linear.blue * fx;
    let bottom_a = c00_linear.alpha * (1.0 - fx) + c10_linear.alpha * fx;

    let top_r = c01_linear.red * (1.0 - fx) + c11_linear.red * fx;
    let top_g = c01_linear.green * (1.0 - fx) + c11_linear.green * fx;
    let top_b = c01_linear.blue * (1.0 - fx) + c11_linear.blue * fx;
    let top_a = c01_linear.alpha * (1.0 - fx) + c11_linear.alpha * fx;

    let final_r = bottom_r * (1.0 - fy) + top_r * fy;
    let final_g = bottom_g * (1.0 - fy) + top_g * fy;
    let final_b = bottom_b * (1.0 - fy) + top_b * fy;
    let final_a = bottom_a * (1.0 - fy) + top_a * fy;

    Color::linear_rgba(final_r, final_g, final_b, final_a)
}
