use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use bevy::render::texture::{ImageSampler, ImageSamplerDescriptor};
use crate::atmosphere::grid::{AtmosphereGrid, AtmosphereCell};

pub const COLORMAP_SIZE: usize = 256;

// Visualization ranges (can later be made runtime-configurable)
pub const T_MIN_VIS: f32 = 200.0;
pub const T_MAX_VIS: f32 = 350.0;
pub const P_MIN_VIS: f32 = 0.0;
pub const P_MAX_VIS: f32 = 200000.0; // ~2 atm
pub const V_MIN_VIS: f32 = 0.0;
pub const V_MAX_VIS: f32 = 50.0; // m/s

#[derive(Resource)]
pub struct HeatmapTexture(pub Handle<Image>);

#[derive(Resource)]
pub struct HeatmapColormap(pub [u8; COLORMAP_SIZE * 4]);

#[derive(Clone, Copy, Debug)]
pub enum HeatmapField {
    Temperature,
    Pressure,
    VelocityMag,
}

#[derive(Resource)]
pub struct ActiveHeatmapField(pub HeatmapField);

/// Marker component for the heatmap sprite entity
#[derive(Component)]
pub struct HeatmapSprite;

/// Colormap generation - blue → cyan → yellow → red gradient
impl FromWorld for HeatmapColormap {
    fn from_world(_world: &mut World) -> Self {
        let mut data = [0u8; COLORMAP_SIZE * 4];

        for i in 0..COLORMAP_SIZE {
            let t = i as f32 / (COLORMAP_SIZE - 1) as f32;

            // Blue → Cyan → Yellow → Red gradient
            // - Red ramps up linearly: r = t
            // - Green peaks in middle: g = (0.8*t)^1.2 for warmer tones
            // - Blue fades out linearly: b = 1 - t
            let r = (t * 255.0) as u8;
            let g = ((t * 0.8).powf(1.2) * 255.0) as u8;
            let b = ((1.0 - t) * 255.0) as u8;

            let idx = i * 4;
            data[idx + 0] = r;
            data[idx + 1] = g;
            data[idx + 2] = b;
            data[idx + 3] = 255;
        }

        HeatmapColormap(data)
    }
}

/// Normalize temperature to [0, 1] range for colormap lookup
#[inline]
pub fn normalize_temp(t: f32) -> f32 {
    ((t - T_MIN_VIS) / (T_MAX_VIS - T_MIN_VIS)).clamp(0.0, 1.0)
}

/// Normalize pressure to [0, 1] range for colormap lookup
#[inline]
pub fn normalize_pressure(p: f32) -> f32 {
    ((p - P_MIN_VIS) / (P_MAX_VIS - P_MIN_VIS)).clamp(0.0, 1.0)
}

/// Normalize velocity magnitude to [0, 1] range for colormap lookup
#[inline]
pub fn normalize_velocity(v: f32) -> f32 {
    ((v - V_MIN_VIS) / (V_MAX_VIS - V_MIN_VIS)).clamp(0.0, 1.0)
}

/// Normalize scalar field value based on field type
#[inline]
pub fn normalize_scalar(scalar: f32, field: HeatmapField) -> f32 {
    match field {
        HeatmapField::Temperature => normalize_temp(scalar),
        HeatmapField::Pressure => normalize_pressure(scalar),
        HeatmapField::VelocityMag => normalize_velocity(scalar),
    }
}

/// Sample scalar field from cell
#[inline]
fn sample_scalar(cell: &AtmosphereCell, field: HeatmapField) -> f32 {
    match field {
        HeatmapField::Temperature => cell.temperature,
        HeatmapField::Pressure    => cell.pressure,
        HeatmapField::VelocityMag => (cell.u * cell.u + cell.v * cell.v).sqrt(),
    }
}

/// Plugin for scalar field heatmap visualization
pub struct HeatmapPlugin;

impl Plugin for HeatmapPlugin {
    fn build(&self, app: &mut App) {
        app
            .init_resource::<HeatmapColormap>()
            .insert_resource(ActiveHeatmapField(HeatmapField::Temperature))
            .add_systems(Startup, setup_heatmap_texture)
            .add_systems(Update, (toggle_heatmap_visibility, update_heatmap_system));
    }
}

/// Create heatmap texture and spawn rendering quad
fn setup_heatmap_texture(
    mut commands: Commands,
    mut images: ResMut<Assets<Image>>,
    atmosphere: Res<AtmosphereGrid>,
) {
    let grid_w = atmosphere.width as u32;
    let grid_h = atmosphere.height as u32;

    let size = Extent3d {
        width: grid_w,
        height: grid_h,
        depth_or_array_layers: 1,
    };

    let mut image = Image::new_fill(
        size,
        TextureDimension::D2,
        &[0, 0, 0, 255],
        TextureFormat::Rgba8UnormSrgb,
        bevy::render::render_asset::RenderAssetUsages::RENDER_WORLD | bevy::render::render_asset::RenderAssetUsages::MAIN_WORLD,
    );

    // Linear filtering for smooth interpolation between cells
    image.sampler = ImageSampler::Descriptor(ImageSamplerDescriptor::linear());

    let handle = images.add(image);
    commands.insert_resource(HeatmapTexture(handle.clone()));

    let world_w = atmosphere.width as f32 * atmosphere.tile_size_physical;
    let world_h = atmosphere.height as f32 * atmosphere.tile_size_physical;

    // Spawn heatmap sprite as a 1x1 normalized quad scaled to world dimensions
    // Render at z=0 (background layer) - tiles/walls should use z > 0 to appear above
    commands.spawn((
        SpriteBundle {
            texture: handle,
            transform: Transform::from_xyz(world_w * 0.5, world_h * 0.5, 0.0)
                .with_scale(Vec3::new(world_w, world_h, 1.0)),
            ..Default::default()
        },
        HeatmapSprite,
    ));

    info!("Heatmap visualization initialized: {}x{} texture", grid_w, grid_h);
}

/// Toggle heatmap visibility with '1' key
fn toggle_heatmap_visibility(
    keyboard: Res<ButtonInput<KeyCode>>,
    mut query: Query<&mut Visibility, With<HeatmapSprite>>,
) {
    if keyboard.just_pressed(KeyCode::Digit1) {
        for mut visibility in query.iter_mut() {
            *visibility = match *visibility {
                Visibility::Hidden => {
                    info!("Texture-based heatmap enabled");
                    Visibility::Visible
                }
                _ => {
                    info!("Texture-based heatmap disabled");
                    Visibility::Hidden
                }
            };
        }
    }
}

/// Update heatmap texture from atmosphere grid each frame
fn update_heatmap_system(
    atmosphere: Res<AtmosphereGrid>,
    mut images: ResMut<Assets<Image>>,
    heatmap: Res<HeatmapTexture>,
    colormap: Res<HeatmapColormap>,
    field: Res<ActiveHeatmapField>,
) {
    let grid_w = atmosphere.width as usize;
    let grid_h = atmosphere.height as usize;

    let image = match images.get_mut(&heatmap.0) {
        Some(img) => img,
        None => return,
    };

    debug_assert_eq!(image.width() as usize, grid_w);
    debug_assert_eq!(image.height() as usize, grid_h);

    let data = &mut image.data;
    let lut = &colormap.0;
    let field = field.0;

    for y in 0..grid_h {
        for x in 0..grid_w {
            let idx_cell = atmosphere.index(x as u32, y as u32);
            let cell = &atmosphere.cells[idx_cell];

            // Sample the scalar field and normalize based on field type
            let scalar = sample_scalar(cell, field);
            let t_norm = normalize_scalar(scalar, field);
            let lut_idx = (t_norm * ((COLORMAP_SIZE - 1) as f32)).round() as usize;

            debug_assert!(lut_idx < COLORMAP_SIZE, "LUT index out of bounds: {}", lut_idx);

            let src = &lut[lut_idx * 4..lut_idx * 4 + 4];
            let idx_px = (y * grid_w + x) * 4;

            data[idx_px + 0] = src[0];
            data[idx_px + 1] = src[1];
            data[idx_px + 2] = src[2];
            data[idx_px + 3] = 255;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_temperature_normalization_clamps() {
        // Below min
        assert_eq!(normalize_temp(150.0), 0.0);
        // Above max
        assert_eq!(normalize_temp(400.0), 1.0);
        // At min
        assert_eq!(normalize_temp(T_MIN_VIS), 0.0);
        // At max
        assert_eq!(normalize_temp(T_MAX_VIS), 1.0);
        // Midpoint
        let mid = (T_MIN_VIS + T_MAX_VIS) / 2.0;
        assert!((normalize_temp(mid) - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_pressure_normalization_clamps() {
        // Below min
        assert_eq!(normalize_pressure(-1000.0), 0.0);
        // Above max
        assert_eq!(normalize_pressure(300000.0), 1.0);
        // At min
        assert_eq!(normalize_pressure(P_MIN_VIS), 0.0);
        // At max
        assert_eq!(normalize_pressure(P_MAX_VIS), 1.0);
        // Midpoint
        let mid = (P_MIN_VIS + P_MAX_VIS) / 2.0;
        assert!((normalize_pressure(mid) - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_velocity_normalization_clamps() {
        // Below min
        assert_eq!(normalize_velocity(-10.0), 0.0);
        // Above max
        assert_eq!(normalize_velocity(100.0), 1.0);
        // At min
        assert_eq!(normalize_velocity(V_MIN_VIS), 0.0);
        // At max
        assert_eq!(normalize_velocity(V_MAX_VIS), 1.0);
        // Midpoint
        let mid = (V_MIN_VIS + V_MAX_VIS) / 2.0;
        assert!((normalize_velocity(mid) - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_normalize_scalar_delegates_correctly() {
        // Temperature
        let t = 275.0;
        assert_eq!(normalize_scalar(t, HeatmapField::Temperature), normalize_temp(t));

        // Pressure
        let p = 100000.0;
        assert_eq!(normalize_scalar(p, HeatmapField::Pressure), normalize_pressure(p));

        // Velocity
        let v = 25.0;
        assert_eq!(normalize_scalar(v, HeatmapField::VelocityMag), normalize_velocity(v));
    }

    #[test]
    fn test_sample_scalar_extracts_correct_field() {
        let cell = AtmosphereCell {
            rho_o2: 0.3,
            rho_n2: 0.9,
            rho_co2: 0.001,
            u: 3.0,
            v: 4.0,
            temperature: 293.0,
            pressure: 101325.0,
        };

        assert_eq!(sample_scalar(&cell, HeatmapField::Temperature), 293.0);
        assert_eq!(sample_scalar(&cell, HeatmapField::Pressure), 101325.0);

        // Velocity magnitude should be sqrt(3^2 + 4^2) = 5.0
        let vel_mag = sample_scalar(&cell, HeatmapField::VelocityMag);
        assert!((vel_mag - 5.0).abs() < 0.001);
    }

    #[test]
    fn test_colormap_endpoints() {
        let mut world = World::new();
        let map = HeatmapColormap::from_world(&mut world);

        // Blue at t=0 (cold)
        let idx_start = 0;
        assert_eq!(map.0[idx_start + 0], 0);   // R = 0
        assert_eq!(map.0[idx_start + 3], 255); // Alpha = 255
        assert!(map.0[idx_start + 2] > 200);   // B should be high

        // Red at t=1 (hot)
        let idx_end = (COLORMAP_SIZE - 1) * 4;
        assert!(map.0[idx_end + 0] > 200); // R should be high
        assert_eq!(map.0[idx_end + 2], 0); // B = 0
        assert_eq!(map.0[idx_end + 3], 255); // Alpha = 255
    }

    #[test]
    fn test_colormap_gradient_monotonic() {
        let mut world = World::new();
        let map = HeatmapColormap::from_world(&mut world);

        // Red channel should increase monotonically
        for i in 0..(COLORMAP_SIZE - 1) {
            let r1 = map.0[i * 4 + 0];
            let r2 = map.0[(i + 1) * 4 + 0];
            assert!(r2 >= r1, "Red should increase: {} -> {}", r1, r2);
        }

        // Blue channel should decrease monotonically
        for i in 0..(COLORMAP_SIZE - 1) {
            let b1 = map.0[i * 4 + 2];
            let b2 = map.0[(i + 1) * 4 + 2];
            assert!(b2 <= b1, "Blue should decrease: {} -> {}", b1, b2);
        }
    }
}
