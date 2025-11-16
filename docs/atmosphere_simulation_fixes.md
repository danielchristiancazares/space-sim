# Atmosphere Simulation Fixes

## Critical Fixes

### 1. Fix Vacuum Handling in Pressure Correction (`src/atmosphere/steps/pressure_projection.rs`)

**Problem:** Clamping density to `1e-6` in vacuum causes incorrect velocity corrections and divergence issues.

**Solution:** Skip pressure correction for vacuum cells:

- Check if `total_density() < MIN_DENSITY_THRESHOLD` before processing
- Set pressure correction to 0 for vacuum cells
- Skip velocity updates for vacuum cells

**Files:** `src/atmosphere/steps/pressure_projection.rs` (lines 61-112, 114-155)

### 2. Improve Boundary Handling in Advection (`src/atmosphere/steps/advection.rs`)

**Problem:** Clamping coordinates to grid bounds causes mass accumulation at walls.

**Solution:** Implement proper boundary conditions:

- For blocked cells: use no-slip boundary (zero velocity)
- For out-of-bounds: use reflective or extrapolation boundary
- Consider adding boundary condition enum for future flexibility

**Files:** `src/atmosphere/steps/advection.rs` (lines 108-149, specifically `interpolate_cell`)

## Medium Priority Fixes

### 3. Extract Magic Numbers to Constants (`src/atmosphere/constants.rs`)

**Problem:** Magic numbers scattered throughout code make tuning difficult.

**Solution:** Add constants:

- `CMB_TEMPERATURE: f32 = 2.7`
- `MAX_TEMPERATURE: f32 = 10000.0`
- `MIN_DENSITY_THRESHOLD: f32 = 1e-6`
- Replace all occurrences in codebase

**Files:**

- `src/atmosphere/constants.rs` (add new constants)
- `src/atmosphere/steps/advection.rs` (lines 75, 98, 101)
- `src/atmosphere/steps/pressure_projection.rs` (line 71)
- `src/atmosphere/grid.rs` (line 75)

### 4. Fix Temperature Handling in Life Support (`src/atmosphere/sources.rs`)

**Problem:** Life support only clamps temperature, doesn't add thermal energy.

**Solution:** Add thermal energy based on gas generation:

- Calculate thermal energy from gas generation rate
- Update temperature using heat capacity
- Keep clamping as safety fallback

**Files:** `src/atmosphere/sources.rs` (lines 36-38)

### 5. Extract Divergence Calculation (`src/atmosphere/steps/pressure_projection.rs`, `src/atmosphere/monitoring.rs`)

**Problem:** Divergence calculation duplicated in two files.

**Solution:** Create shared helper function:

- Extract to `grid.rs` or new `utils.rs`
- Use in both `pressure_projection.rs` and `monitoring.rs`
- Ensure consistent boundary handling

**Files:**

- `src/atmosphere/steps/pressure_projection.rs` (lines 22-58)
- `src/atmosphere/monitoring.rs` (lines 83-122)
- Location for helper function (new or existing module)

## Low Priority Improvements

### 6. Add Validation Helpers (`src/atmosphere/grid.rs`)

**Problem:** No validation for invalid cell states.

**Solution:** Add `is_valid()` method to `AtmosphereCell`:

- Check for negative densities
- Check for NaN/Inf in velocities
- Check for valid temperature range
- Optionally add `validate()` method returning errors

**Files:** `src/atmosphere/grid.rs` (extend `AtmosphereCell` impl)

### 7. Document Pressure Update Sequence (`src/atmosphere/simulation.rs`)

**Problem:** The paired `update_grid_pressures` calls appear redundant, but both are required so that `pressure_flux_step` and following systems operate on pressures that match the latest densities. Without documentation, this looks like unnecessary work and invites regressions.

**Solution:**

- Keep both `update_grid_pressures` calls.
- Add in-line comments explaining why each call is required.
- If future optimization is attempted, ensure the flux stage still sees freshly computed pressures and that cached pressures remain synchronized afterward.

**Files:** `src/atmosphere/simulation.rs` (lines 34-45)

### 8. Improve Diffusion Boundary Conditions (`src/atmosphere/steps/diffusion.rs`)

**Problem:** Using center cell values at boundaries may not match intended physics.

**Solution:**

- Document the current boundary choice
- Optionally switch to explicit no-flux (zero gradient) boundary handling

**Files:** `src/atmosphere/steps/diffusion.rs` (lines 96-131)

## Testing

### 9. Add Tests for Critical Fixes

- Mass conservation with pressure flux + vacuum cells
- Advection boundary conditions and vacuum expansion
- Temperature updates from life support

**Files:** Extend existing test modules or add new ones under `src/atmosphere/steps/tests` and `src/atmosphere/sources/tests`

