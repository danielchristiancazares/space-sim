use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use bevy::render::texture::{ImageSampler, ImageSamplerDescriptor};
use crate::atmosphere::grid::{AtmosphereGrid, AtmosphereCell};

pub const COLORMAP_SIZE: usize = 256;

// Visualization temperature range (can later be made runtime-configurable)
pub const T_MIN_VIS: f32 = 200.0;
pub const T_MAX_VIS: f32 = 350.0;

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

/// Colormap generation - blue → cyan → yellow → red gradient
impl FromWorld for HeatmapColormap {
    fn from_world(_world: &mut World) -> Self {
        let mut data = [0u8; COLORMAP_SIZE * 4];

        for i in 0..COLORMAP_SIZE {
            let t = i as f32 / (COLORMAP_SIZE - 1) as f32;

            // Blue → cyan → yellow → red gradient
            let r = (t.powf(1.0) * 255.0) as u8;
            let g = ((t * 0.8).powf(1.2) * 255.0) as u8;
            let b = ((1.0 - t).powf(1.0) * 255.0) as u8;

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

/// Sample scalar field from cell (for future multi-field support)
#[inline]
#[allow(dead_code)]
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
            .add_systems(Update, update_heatmap_system);
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

    // Center the quad over the simulation domain
    commands.spawn(SpriteBundle {
        texture: handle,
        transform: Transform::from_xyz(world_w * 0.5, world_h * 0.5, 0.0)
            .with_scale(Vec3::new(world_w, world_h, 1.0)),
        ..Default::default()
    });

    info!("Heatmap visualization initialized: {}x{} texture", grid_w, grid_h);
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

            let scalar = match field {
                HeatmapField::Temperature => cell.temperature,
                HeatmapField::Pressure    => cell.pressure,
                HeatmapField::VelocityMag => (cell.u * cell.u + cell.v * cell.v).sqrt(),
            };

            // For now reuse temperature normalization; for other fields
            // add dedicated normalization functions.
            let t_norm = normalize_temp(scalar);
            let lut_idx = (t_norm * ((COLORMAP_SIZE - 1) as f32)).round() as usize;

            let src = &lut[lut_idx * 4..lut_idx * 4 + 4];
            let idx_px = (y * grid_w + x) * 4;

            data[idx_px + 0] = src[0];
            data[idx_px + 1] = src[1];
            data[idx_px + 2] = src[2];
            data[idx_px + 3] = 255;
        }
    }
}
