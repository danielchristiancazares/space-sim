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

#[derive(Clone, Copy, Debug)]
pub struct VisualizationState {
    heatmap_enabled: bool,
    arrows_enabled: bool,
}

impl Default for VisualizationState {
    fn default() -> Self {
        Self {
            heatmap_enabled: true,
            arrows_enabled: false,
        }
    }
}

impl VisualizationState {
    fn reset(&mut self) {
        *self = Self::default();
    }
}

pub fn debug_visualization(
    atmosphere: Res<AtmosphereGrid>,
    collision_map: Res<TileCollisionMap>,
    keyboard: Res<ButtonInput<KeyCode>>,
    mut state: Local<VisualizationState>,
    mut gizmos: Gizmos,
) {
    if keyboard.just_pressed(KeyCode::Digit3) {
        state.reset();
        info!("Atmosphere visualization reset (velocity arrows off)");
        info!("Note: Press '1' to toggle texture-based heatmap (handled by HeatmapPlugin)");
    }

    // Note: '1' key toggle is now handled by HeatmapPlugin (texture-based heatmap)
    // This gizmo-based visualization is deprecated
    if keyboard.just_pressed(KeyCode::Digit1) {
        state.heatmap_enabled = !state.heatmap_enabled;
        // No longer log here - HeatmapPlugin handles the toggle
    }

    if keyboard.just_pressed(KeyCode::Digit2) {
        state.arrows_enabled = !state.arrows_enabled;
        info!(
            "Atmosphere velocity arrows {}",
            if state.arrows_enabled {
                "enabled"
            } else {
                "disabled"
            }
        );
    }

    if !state.heatmap_enabled && !state.arrows_enabled {
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

            // Note: Gizmo-based heatmap is deprecated in favor of texture-based heatmap
            // This code is kept for reference but not rendered
            // The texture-based heatmap (HeatmapPlugin) now handles visualization
            if false && state.heatmap_enabled {
                let color = get_cell_color(
                    &atmosphere,
                    &collision_map,
                    x,
                    y,
                    VisualizationMode::Breathability,
                );
                gizmos.rect_2d(
                    Isometry2d::from_translation(Vec2::new(world_x, world_y)),
                    Vec2::splat(tile_size * 0.95),
                    color,
                );
            }

            if state.arrows_enabled {
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
