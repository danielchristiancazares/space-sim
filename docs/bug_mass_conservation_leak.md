# Bug Report: Mass Conservation Violation in Atmospheric Simulation

**Date**: 2025-11-12
**Severity**: Critical
**Status**: Partially Fixed (88.5% mass loss remaining)
**Test**: `life_support_fill_remains_stable` (tests/atmosphere_integration.rs)

---

## Summary

The atmospheric simulation loses 88.5% of total mass over a 15-minute simulation period. Life support systems inject 121.8 kg of gas, but only 14.1 kg remains at completion. This violates fundamental conservation laws and makes the simulation unusable for realistic gameplay.

---

## Test Results

### Before Fixes (2025-11-12 22:00 PST)

```
Test: life_support_fill_remains_stable
Status: FAILED (numerical explosion)
Runtime: 403 seconds (6.7 minutes) before crash
Error: pressure_correction_step panicked
  rhs = -2.05×10³⁷ Pa (catastrophic overflow)
  neighbors = (8.13×10³⁷, 8.36×10³⁷, ...)
Final mass: 0.785 kg / 121.800 kg expected
Mass loss: 99.4%
```

### After Partial Fixes (2025-11-12 01:30 PST)

```
Test: life_support_fill_remains_stable
Status: FAILED (mass conservation violation)
Runtime: 1100 seconds (18.3 minutes) - completed without crash ✓
Error: assertion failed on final mass check
Final mass: 14.074 kg / 121.800 kg expected (±5% tolerance)
Mass loss: 88.5%
Improvement: 18× reduction in mass loss (0.785 → 14.074 kg)
```

---

## Root Cause Analysis

### Primary Issue: Mass Leakage Through Multiple Mechanisms

The simulation pipeline has three stages that should conserve mass:

1. **Advection** (MacCormack) - ✅ Conservative by design
2. **Diffusion** (Laplacian) - ❌ Was leaking at boundaries (FIXED)
3. **Pressure Flux** - ⚠️ Suspected leak in negative density clamping
4. **Pressure Projection** - ✅ Only modifies velocity (not mass)

### Issues Fixed

#### 1. Vacuum Threshold Too Small (P0 - FIXED)

**Location**: `src/atmosphere/steps/pressure_flux.rs:52,62,94,128`

**Problem**:
```rust
// BEFORE: Using f32::EPSILON ≈ 1.2×10⁻⁷
if donor_state.pressure <= f32::EPSILON {
```

For pressure in Pascals:
- Earth atmosphere: ~10⁵ Pa
- Hard vacuum: ~10⁻⁶ Pa
- f32::EPSILON: ~10⁻⁷ (meaningless for pressure)

Cells with 1 Pa pressure (practically vacuum) were treated as normal, causing:
- Division by near-zero in momentum calculations
- Incorrect flux calculations
- Potential NaN propagation

**Fix**:
```rust
// AFTER: Physically meaningful threshold
const VACUUM_THRESHOLD_PA: f32 = 1.0;
if donor_state.pressure <= VACUUM_THRESHOLD_PA {
```

**Impact**: Prevents numerical instability in vacuum regions.

---

#### 2. Incorrect Wall Boundary Conditions in Divergence (P1 - FIXED)

**Location**: `src/atmosphere/steps/pressure_projection.rs:30-56`

**Problem**:
```rust
// BEFORE: Used center cell velocity for walls
let u_right = if wall_or_boundary {
    atmosphere.cells[idx].u  // Wrong: assumes u_wall = u_center
} else {
    neighbor.u
};
```

This treats walls as frictionless (free-slip) instead of no-slip.

**Correct physics**:
- At walls: velocity = 0 (no-slip boundary condition)
- Divergence should account for wall at u=0, not center velocity

**Fix**:
```rust
// AFTER: No-slip boundary condition
let u_right = if wall_or_boundary {
    0.0  // Wall velocity (no-slip)
} else {
    neighbor.u
};
```

**Impact**: Correctly computes divergence near walls, preventing Poisson solver divergence.

---

#### 3. Diffusion Boundary Treated as Vacuum (P0 - FIXED)

**Location**: `src/atmosphere/steps/diffusion.rs:112-122`

**Problem**:
```rust
// BEFORE: Grid edges treated as vacuum
} else {
    AtmosphereCell {
        rho_o2: 0.0,   // Vacuum
        rho_n2: 0.0,
        rho_co2: 0.0,
        ...
    }
}
```

When computing Laplacian at edge cells:
```
∇²ρ = (ρ_L + ρ_R + ρ_U + ρ_D - 4ρ_C) / dx²
```

If boundary returns 0 (vacuum), Laplacian becomes negative:
```
∇²ρ ≈ 0 + neighbor + neighbor + 0 - 4×ρ_C = negative
ρ_new = ρ_old + D×dt×∇²ρ = mass diffuses out
```

**Fix**:
```rust
// AFTER: Zero-gradient boundary (closed room)
} else {
    AtmosphereCell {
        rho_o2: center.rho_o2,  // Same as center (no gradient)
        rho_n2: center.rho_n2,
        rho_co2: center.rho_co2,
        ...
    }
}
```

Now Laplacian at edges = 0 (no flux across boundary).

**Impact**: Eliminates 10% of mass loss, but 88% still remains.

---

### Issues Remaining (Suspected)

#### 4. Negative Density Clamping (P0 - NEEDS INVESTIGATION)

**Location**: `src/atmosphere/steps/pressure_flux.rs:183-187`

```rust
// Non-conservative safety net
for cell in atmosphere.cells.iter_mut() {
    cell.rho_o2 = cell.rho_o2.max(0.0);
    cell.rho_n2 = cell.rho_n2.max(0.0);
    cell.rho_co2 = cell.rho_co2.max(0.0);
}
```

**Why this exists**: The availability scaling (lines 158-177) should prevent negative densities. If densities go negative despite this, it indicates a bug in the flux calculation.

**Why it's bad**: Clamping negative values to zero **destroys mass**. If a cell has -5 kg/m³ O2 (unphysical), clamping to 0 deletes that mass from the universe.

**Hypothesis**: This is where the remaining 88% is being lost.

**Evidence needed**:
- Log total negative mass before clamping
- Track spatial distribution of negative cells
- Track how often this safety net triggers
- Determine if it's numerical precision or logic error

**Leading hypothesis**: Availability scaling has subtle bug allowing slightly more flux than available. Over 18 minutes (thousands of timesteps), small errors compound into massive negative densities that get clamped away, destroying 88% of total mass.

---

#### 5. Advection Boundary Handling (P2 - SUSPECTED)

**Location**: `src/atmosphere/steps/advection.rs:92-98`

```rust
let x0 = x.floor().max(0.0).min((width - 1) as f32) as u32;
let y0 = y.floor().max(0.0).min((height - 1) as f32) as u32;
```

When backtracing goes outside grid:
- Clamps to edge cells
- Doesn't distinguish walls vs open boundaries
- May interpolate from vacuum incorrectly

**Impact**: Minor (most particles stay in grid due to CFL limiting), but could contribute to leak.

---

#### 6. Pressure Correction at Vacuum Interface (P2 - SUSPECTED)

**Location**: `src/atmosphere/steps/pressure_projection.rs:129-165`

```rust
if rho < 1e-6 {
    continue;  // Skip vacuum cells
}
```

Velocity corrections skip vacuum cells, but what about cells adjacent to vacuum?

**Possible issue**: Momentum may not be conserved when correcting velocities near vacuum boundaries.

---

## Reproduction Steps

```bash
cargo test --test atmosphere_integration life_support_fill_remains_stable
```

**Expected**: Final mass 121.8 kg ± 5%
**Actual**: Final mass 14.1 kg (88.5% loss)
**Runtime**: ~18 minutes

---

## Diagnostic Recommendations

### 1. Add Mass Leak Tracking with Spatial Distribution (High Priority)

Add logging before negative clamping in `pressure_flux_step`:

```rust
// BEFORE clamping - Track spatial distribution
let mut negative_cells = 0;
let mut total_neg_mass = 0.0;
let mut max_negative = 0.0f32;

for (idx, cell) in atmosphere.cells.iter().enumerate() {
    let neg_o2 = cell.rho_o2.min(0.0).abs();
    let neg_n2 = cell.rho_n2.min(0.0).abs();
    let neg_co2 = cell.rho_co2.min(0.0).abs();
    let neg_total = neg_o2 + neg_n2 + neg_co2;

    if neg_total > 0.0 {
        negative_cells += 1;
        total_neg_mass += neg_total;
        max_negative = max_negative.max(neg_total);

        // Log first 10 negative cells to identify spatial pattern
        if negative_cells <= 10 {
            let x = idx % atmosphere.width as usize;
            let y = idx / atmosphere.width as usize;
            eprintln!(
                "Negative density at ({}, {}): O2={:.3e}, N2={:.3e}, CO2={:.3e}",
                x, y, cell.rho_o2, cell.rho_n2, cell.rho_co2
            );
        }
    }
}

if total_neg_mass > 1e-6 {
    eprintln!(
        "WARNING: Clamping {:.3e} kg of negative mass across {} cells (max: {:.3e} kg/m³)",
        total_neg_mass, negative_cells, max_negative
    );
}

// AFTER clamping
for cell in atmosphere.cells.iter_mut() {
    cell.rho_o2 = cell.rho_o2.max(0.0);
    cell.rho_n2 = cell.rho_n2.max(0.0);
    cell.rho_co2 = cell.rho_co2.max(0.0);
}
```

**Why spatial distribution matters**:
- **Clustered negatives** → Interface problem (vacuum boundaries, flux discontinuities)
- **Scattered negatives** → Numerical precision issue (accumulation over many timesteps)
- **Edge-only negatives** → Boundary condition bug
- **Random distribution** → Availability scaling bug allowing over-extraction

### 2. Track Mass Through Pipeline (High Priority)

In `simulate_atmosphere` (simulation.rs), log mass before/after each step:

```rust
fn simulate_atmosphere(...) {
    let mass_before = compute_total_mass(&atmosphere);

    advection_step(&mut atmosphere, &collision_map, dt_sim);
    let mass_after_advection = compute_total_mass(&atmosphere);

    diffusion_step(&mut atmosphere, &collision_map, dt_sim);
    let mass_after_diffusion = compute_total_mass(&atmosphere);

    // ... etc

    eprintln!("Mass: before={}, advection={}, diffusion={}, flux={}, projection={}",
        mass_before, mass_after_advection, mass_after_diffusion, ...);
}
```

### 3. Check Availability Scaling (Medium Priority)

Add assertion in `process_flux` before mass transfer:

```rust
let total_available = avail_o2 + avail_n2 + avail_co2;
let total_requested = m_o2 + m_n2 + m_co2;
assert!(
    total_requested <= total_available * 1.001,
    "Requesting {} kg but only {} kg available",
    total_requested, total_available
);
```

---

## Timeline

| Time | Event |
|------|-------|
| 22:00 PST | Initial test run: 99.4% mass loss, crash at 403s |
| 22:30 PST | Code review identified vacuum threshold bug |
| 01:00 PST | Applied fixes for vacuum threshold, wall divergence, diffusion boundaries |
| 01:30 PST | Retest: 88.5% mass loss, stable run 1100s |

---

## Related Files

- `src/atmosphere/steps/pressure_flux.rs:183-187` - Negative clamping (suspected leak)
- `src/atmosphere/steps/diffusion.rs:112-123` - Boundary conditions (fixed)
- `src/atmosphere/steps/pressure_projection.rs:30-56` - Divergence calc (fixed)
- `src/atmosphere/steps/advection.rs:92-98` - Boundary interpolation (suspected)
- `src/atmosphere/simulation.rs:12-60` - Main loop (add diagnostics here)

---

## Next Steps

1. **Add mass leak diagnostics** to identify which step loses mass
2. **Investigate negative density clamping** - log how much mass is destroyed
3. **Check availability scaling** - verify flux calculations are correct
4. **Review advection boundaries** - ensure backtracing doesn't lose mass
5. **Add per-step mass conservation tests** - isolate the leak

---

## Expected Behavior

For a **closed room** with life support:
- Mass injected by life support = mass retained (±1% numerical error)
- No mass should escape through walls or grid boundaries
- Conservation error should be < 1% over 15 minutes

For **open to space**:
- Mass can escape through designated openings (airlocks, breaches)
- Leak rate should match physical flow rate
- Currently all boundaries are closed, so 0 mass should escape

---

## Workarounds

None. Mass conservation is fundamental. Simulation is currently non-physical.

---

## Additional Notes

- Unit tests (7/7) all pass - conservation works for simple cases
- Integration tests show issue only appears over long timescales
- Gradual accumulation suggests numerical drift, not discrete bug
- 18-minute test runtime makes iteration slow

**Recommendation**: Add extensive logging, then binary search to isolate which step(s) lose mass.
