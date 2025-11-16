# Atmospheric Simulation Pipeline

This directory contains the numerical methods that evolve the atmospheric state over time. Each "step" is a substep in the operator splitting approach used to integrate the full Navier-Stokes-like equations.

## Quick Reference

```
Main Loop (../simulation.rs) → Subdivided by CFL condition
    ↓
[1] Advection            → Flux-conservative densities + MacCormack (u, v, T)
[2] Diffusion            → Momentum/thermal/species smoothing
[3] Compression Heating  → Apply -(γ-1)T∇·u
[4] Update Pressures     → Ideal gas law for diagnostics & sources
```

Legacy notes: Earlier experiments added pressure-driven flux and a Poisson-based projection. Those systems are currently disabled and described later in this document for historical context only.

## Detailed Pipeline Flow

```
Main Loop (../simulation.rs:simulate_atmosphere)
    ↓
┌──────────────────────────────────────────────────────┐
│  Per-Frame Loop (subdivided by CFL condition)       │
│  dt_sub = min(remaining, CFL × dx / u_max)          │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [1] ADVECTION (advection.rs)                        │
├──────────────────────────────────────────────────────┤
│  Method: MacCormack (2nd-order)                     │
│  ├─ Forward:  φ* = φⁿ backtraced by u×dt            │
│  ├─ Backward: φ** = φ* forwardtraced by u×dt        │
│  └─ Correct:  φⁿ⁺¹ = φ* + 0.5(φⁿ - φ**)             │
│  Updates: rho_o2, rho_n2, rho_co2, u, v, T          │
│  Buffers: cells_buffer → forward → backward → cells │
│  Purpose: Transport quantities along flow           │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [2] DIFFUSION (diffusion.rs)                        │
├──────────────────────────────────────────────────────┤
│  Method: Explicit Euler on Laplacian                │
│  ∂φ/∂t = α∇²φ                                        │
│  ∇²φ = (L + R + U + D - 4C) / dx²                   │
│  Three processes:                                    │
│  ├─ Momentum:     ν = 1.5×10⁻⁵ m²/s (viscosity)     │
│  ├─ Thermal:      α = 2.1×10⁻⁵ m²/s (heat)          │
│  └─ Species:      D = 5.0×10⁻⁴ m²/s (diffusion)     │
│  Stability: α×dt/dx² ≤ 0.2 (enforced)               │
│  Purpose: Smooth sharp gradients                    │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [3] UPDATE PRESSURES (simulation.rs)                │
├──────────────────────────────────────────────────────┤
│  Ideal gas: P = (ρ/M) × R × T                       │
│  For each species: P_i = (ρ_i/M_i) × R × T          │
│  Total: P = P_O2 + P_N2 + P_CO2                      │
│  Purpose: Sync pressure with density/temperature    │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [4] PRESSURE FLUX (pressure_flux.rs) ⚠️ KEY STEP    │
├──────────────────────────────────────────────────────┤
│  Method: Donor cell with momentum conservation      │
│  For each cell pair:                                 │
│  1. Compute ΔP = P_high - P_low                      │
│  2. Calculate m_flux = k × ΔP × area × dt            │
│  3. Determine donor (high P) & receiver (low P)      │
│  4. Split by species: m_i = m × (P_i / P_total)     │
│  5. Check availability: scale if insufficient        │
│  6. Limit: cap at 20% of donor mass                  │
│  7. Transfer mass: donor.ρ -= Δρ, receiver.ρ += Δρ  │
│  8. Transfer momentum: conserve p = m×v              │
│  Updates: rho_o2, rho_n2, rho_co2, u, v             │
│  Purpose: Compressible flow (density changes)       │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [5] UPDATE PRESSURES (simulation.rs)                │
├──────────────────────────────────────────────────────┤
│  Recalculate P after mass transfer                  │
│  Prepares state for pressure projection              │
└──────────────────────────────────────────────────────┘
    ↓
┌──────────────────────────────────────────────────────┐
│ [6] PRESSURE PROJECTION (pressure_projection.rs)    │
├──────────────────────────────────────────────────────┤
│  Method: Poisson solve + velocity correction         │
│  Phase 1: Solve ∇²p = (ρ/dt)∇·u                     │
│   └─ Jacobi iteration (max 50, tol 1e-4)            │
│  Phase 2: Correct u' = u - (dt/ρ)∇p                 │
│  Updates: u, v (removes divergence)                  │
│  Purpose: Enforce incompressibility ∇·u = 0          │
└──────────────────────────────────────────────────────┘
    ↓
    Repeat until frame dt consumed
```

---

## Step 1: Advection (advection.rs)

### Purpose
Transports all quantities (densities, velocities, temperature) along the velocity field.

### Algorithm: MacCormack (2nd-order)

The MacCormack scheme is a predictor-corrector method that achieves 2nd-order accuracy:

**Step 1 - Forward (Predictor)**:
```
x_source = x_current - u × dt
φ_forward = interpolate(φⁿ, x_source)
```
Backtrace from current position to find where the fluid came from.

**Step 2 - Backward (Corrector)**:
```
x_forward = x_current + u_forward × dt
φ_backward = interpolate(φ_forward, x_forward)
```
Forward-trace from the predicted state.

**Step 3 - Correction**:
```
φⁿ⁺¹ = φ_forward + 0.5 × (φⁿ - φ_backward)
```
The error `(φⁿ - φ_backward)` estimates truncation error; adding half of it improves accuracy.

### Why MacCormack?

**Numerical diffusion comparison**:
- Semi-Lagrangian (1st-order): Sharp features blur quickly
- MacCormack (2nd-order): Sharp features stay sharper longer
- Essential for maintaining distinct gas pockets and vortices

### Transported Quantities
- `rho_o2`, `rho_n2`, `rho_co2` - Gas densities
- `u`, `v` - Velocity components
- `temperature` - Thermal energy

### Conservation
✅ Conservative by construction - mass/momentum only redistributed, not created/destroyed

### Edge Cases
- Bilinear interpolation at grid boundaries
- Clamps to grid edges (no extrapolation)
- Non-negative density enforcement via `.max(0.0)`

```rust
// Example usage
advection_step(&mut atmosphere, &collision_map, dt);
```

---

## Step 2: Diffusion (diffusion.rs)

### Purpose
Smooths quantities due to molecular diffusion (viscosity, heat conduction, species mixing).

### Algorithm: Explicit Euler on Laplacian

Discretized diffusion equation:
```
∂φ/∂t = α∇²φ

∇²φ ≈ (φ_L + φ_R + φ_U + φ_D - 4φ_C) / dx²

φⁿ⁺¹ = φⁿ + α × dt × ∇²φ
```

Where:
- `φ_L, φ_R, φ_U, φ_D` - Left, right, up, down neighbors
- `φ_C` - Center cell
- `α` - Diffusion coefficient [m²/s]

### Three Diffusion Processes

| Process | Coefficient | Symbol | Physical Meaning |
|---------|-------------|--------|------------------|
| Momentum diffusion | 1.5×10⁻⁵ m²/s | ν | Kinematic viscosity (smooths velocity) |
| Thermal diffusion | 2.1×10⁻⁵ m²/s | α_th | Thermal diffusivity (smooths temperature) |
| Species diffusion | 5.0×10⁻⁴ m²/s | D | Gas diffusivity (smooths concentrations) |

### Stability Condition

For explicit 2D diffusion:
```
α × dt / dx² ≤ 0.25  (theoretical limit)
```

Code enforces `MAX_ALPHA = 0.2` (safety margin).

### Boundary Conditions

**Walls** (solid boundaries):
- Velocity: No-slip (u=0, v=0)
- Scalars: Zero-gradient (matches center cell)

**Grid edges** (open boundaries):
- Currently: Vacuum (ρ=0, P=0) → ⚠️ causes mass loss
- Fix: Use zero-gradient for closed rooms

```rust
// Example usage
diffusion_step(&mut atmosphere, &collision_map, dt);
```

---

## Step 3: Compression Heating (compression.rs)

### Purpose
Applies the thermodynamic coupling between velocity divergence and temperature, ensuring compression heats the gas while expansion cools it.

### Governing Equation
```
∂T/∂t = -(γ - 1) T ∇·u
```
with γ sourced from `constants::GAMMA`.

### Implementation Highlights
- Divergence uses central differences with collision-aware sampling (no-slip boundaries around blocked tiles).
- Backward-Euler style update `T_new = T_old / (1 + (γ-1) div_u dt)` prevents runaway amplification under strong divergence.
- Temperatures are clamped between 2.7 K and 10 000 K for stability, and extreme divergence values are capped before applying the update.

### Practical Impact
- Life-support jets or vents that compress a region locally will heat it without relying on additional source systems.
- When air expands into vacuum (positive divergence) the step cools it, matching intuitive “draft feels cold” behavior.

```rust
compression_heating_step(&mut atmosphere, &collision_map, dt);
```

---

## Legacy: Pressure Flux (pressure_flux.rs)

> **Status:** Removed from the runtime loop. The notes below remain for anyone revisiting the donor-cell pressure flux solver.

### Purpose
Mass flows from high-pressure to low-pressure cells, carrying momentum and preserving species composition.

### Physical Model

**Flux formula** (Darcy-like):
```
mass_flux = conductance × ΔP × area × dt
```

Where:
- `conductance = 0.01 kg/(Pa·m·s)` - Empirical flow coefficient
- `ΔP = P_high - P_low` [Pa]
- `area = dx` [m] (interface width per unit depth)
- `dt` [s]

### Algorithm (Donor Cell Method)

For each cell pair (horizontal + vertical):

**1. Compute mass transfer**
```rust
let p_diff = P_left - P_right;
let mass = (conductance * p_diff * dx * dt).abs();
```

**2. Determine donor and receiver**
```rust
let (donor, receiver) = if p_diff > 0.0 {
    (left, right)  // Left is high pressure
} else {
    (right, left)  // Right is high pressure
};
```

**3. Split by species fractions**
```rust
// Partial pressures (ideal gas)
P_O2  = (ρ_O2 / M_O2) × R × T
P_N2  = (ρ_N2 / M_N2) × R × T
P_CO2 = (ρ_CO2 / M_CO2) × R × T

// Species fractions
f_O2  = P_O2 / P_total
f_N2  = P_N2 / P_total
f_CO2 = P_CO2 / P_total

// Mass splits
m_O2  = mass × f_O2
m_N2  = mass × f_N2
m_CO2 = mass × f_CO2
```

**4. Check availability** (never take more than exists)
```rust
scale_O2  = min(available_O2 / m_O2, 1.0)
scale_N2  = min(available_N2 / m_N2, 1.0)
scale_CO2 = min(available_CO2 / m_CO2, 1.0)

availability = min(scale_O2, scale_N2, scale_CO2)

m_O2  *= availability
m_N2  *= availability
m_CO2 *= availability
```

**5. Apply stability limit** (prevents complete evacuation)
```rust
max_mass = donor.total_mass() × MAX_FRACTION  // 20%
if transferred_mass > max_mass {
    scale_down(m_O2, m_N2, m_CO2)
}
```

**6. Transfer mass**
```rust
donor.rho_O2    -= m_O2 / volume
donor.rho_N2    -= m_N2 / volume
donor.rho_CO2   -= m_CO2 / volume

receiver.rho_O2  += m_O2 / volume
receiver.rho_N2  += m_N2 / volume
receiver.rho_CO2 += m_CO2 / volume
```

**7. Transfer momentum** (conserve p = m×v)

*Donor velocity update*:
```rust
m_target = m_before - m_transferred
m_actual = donor.total_mass()  // After transfer
v_after  = (m_target × v_before) / m_actual
```

*Receiver velocity update*:
```rust
p_before = m_receiver × v_receiver
p_added  = m_transferred × v_donor
p_after  = p_before + p_added
v_after  = p_after / m_after
```

### Momentum Conservation Diagram

```text
    Donor (before)          Receiver (before)
    m = 10 kg               m = 5 kg
    v = 3 m/s               v = 1 m/s
    p = 30 kg·m/s           p = 5 kg·m/s
           │                      │
           └──── 2 kg @ 3 m/s ────┘
           │                      │
    Donor (after)           Receiver (after)
    m = 8 kg                m = 7 kg
    p = 24 kg·m/s           p = 11 kg·m/s
    v = 3.0 m/s             v = 1.57 m/s

    Total momentum: 35 kg·m/s (conserved ✓)
```

### Conservation Properties
- ✅ **Mass**: Exact (donor loss = receiver gain)
- ✅ **Momentum**: Explicit balance (p = m×v conserved)
- ✅ **Species**: Preserved via partial pressure weighting

### Key Insight
This step implements **compressible flow** - density changes drive pressure changes, which drive more flow. It's philosophically distinct from the pressure projection step (incompressible).

```rust
// Example usage
pressure_flux_step(&mut atmosphere, &collision_map, dt);
```

---

## Legacy: Pressure Projection (pressure_projection.rs)

> **Status:** Also removed from the live pipeline. This section documents the former Poisson-solve approach in case the projection is reintroduced later.

### Purpose
Removes divergence from velocity field to enforce incompressibility constraint `∇·u = 0`.

### Theoretical Background

**Hodge Decomposition**: Any vector field can be split:
```
u = u_incompressible + ∇φ
```
Where `∇·u_incompressible = 0`.

We compute `φ` (pressure correction) and subtract its gradient to get divergence-free velocity.

### Algorithm (Two Phases)

**Phase 1: Solve Poisson Equation**
```
∇²p = (ρ/dt) × ∇·u
```

Discretized:
```
(p_L + p_R + p_U + p_D - 4p_C) / dx² = (ρ/dt) × div

p_C = 0.25 × (p_L + p_R + p_U + p_D - dx² × ρ × div / dt)
```

Solved via **Jacobi iteration**:
- Max iterations: 50
- Tolerance: 1e-4 (pressure change per iteration)
- Boundary: Neumann (∂p/∂n = 0 at walls)

**Phase 2: Velocity Correction**
```
u_new = u_old - (dt/ρ) × ∇p
```

Discretized:
```
u_new = u_old - (dt/ρ) × (p_R - p_L) / (2×dx)
v_new = v_old - (dt/ρ) × (p_U - p_D) / (2×dx)
```

### Vacuum Handling

Cells with ρ < 1e-6 kg/m³ are skipped:
- Poisson solver: Set pressure correction to 0
- Velocity correction: Skip (avoids division by zero)

### Why Needed?

Navier-Stokes equations for **incompressible flow** require `∇·u = 0`. Advection and diffusion can introduce divergence (compression/expansion), so we project back to divergence-free space.

### Tension with Pressure Flux

- **Pressure Flux**: Treats flow as compressible (density varies)
- **Pressure Projection**: Enforces incompressibility (density constant)

This works for **low Mach number flows** (velocities << 340 m/s), where compressibility effects are small but not negligible.

```rust
// Example usage
pressure_correction_step(&mut atmosphere, &collision_map, dt);
```

---

## Operator Splitting

The full system:
```
∂ρ/∂t + ∇·(ρu) = 0                     (continuity)
∂(ρu)/∂t + ∇·(ρu⊗u) = -∇p + μ∇²u       (momentum)
∂(ρE)/∂t + ∇·((ρE+p)u) = ∇·(k∇T)       (energy)
p = ρRT/M                                (state)
```

These are coupled nonlinear PDEs. **Operator splitting** decouples them:

1. Advection: `∂φ/∂t + u·∇φ = 0`
2. Diffusion: `∂φ/∂t = D∇²φ`
3. Pressure: `∂u/∂t = -∇p/ρ`

Each substep is simpler, but we incur **splitting error** (1st-order). For small dt, error is negligible.

---

## Timestep Subdivision (CFL Condition)

Main loop subdivides each frame based on **Courant-Friedrichs-Lewy (CFL) condition**:

```rust
dt_max = CFL_NUMBER × dx / u_max
```

Where:
- `CFL_NUMBER = 0.5` (safety factor)
- `dx` - Grid spacing [m]
- `u_max` - Maximum velocity magnitude [m/s]

**Physical meaning**: Fluid parcels shouldn't travel more than ~0.5 cells per timestep.

**Why?** Explicit schemes become unstable if `u×dt > dx`.

**Adaptive behavior**:
- High velocities → many small substeps
- Low velocities → fewer large substeps (up to frame dt)

---

## Conservation Summary

| Quantity | Advection | Diffusion | Flux | Projection | Overall |
|----------|-----------|-----------|------|------------|---------|
| **Mass** | ✅ Conserved | ⚠️ Leaks | ✅ Conserved | ➖ N/A | ⚠️ Leaks |
| **Momentum** | ✅ Conserved | ✅ Conserved | ✅ Conserved | ⚠️ Modified | ✅ OK |
| **Energy** | ➖ Not tracked | ➖ Not tracked | ➖ Not tracked | ➖ Not tracked | ➖ N/A |

**Mass leak**: Diffusion at grid boundaries treated as vacuum → mass escapes.
**Fix**: Use zero-gradient boundaries for closed rooms.

**Momentum**: Projection modifies momentum to enforce incompressibility (physically correct for incompressible flow).

---

## Coordinate System & Indexing

```
Y (up)
↑
│
└────→ X (right)

Grid layout (3×3 example):
[0,0]───[1,0]───[2,0]
  │       │       │
[0,1]───[1,1]───[2,1]
  │       │       │
[0,2]───[1,2]───[2,2]

Index: idx = y × width + x
```

**Velocity**:
- `u` - X velocity [m/s] (positive = rightward)
- `v` - Y velocity [m/s] (positive = upward)

---

## Common Pitfalls

### 1. Buffer vs Live State
❌ **Wrong**: Read from `cells` during flux computation
✅ **Right**: Read from `cells_buffer`, write to `cells`

Without buffering, order matters (later cells see modified neighbors).

### 2. Vacuum Thresholds
❌ **Wrong**: `if pressure < f32::EPSILON`
✅ **Right**: `if pressure < 1.0` (Pa)

`f32::EPSILON ≈ 1.2×10⁻⁷` is meaningless for pressure in Pascals.

### 3. Pressure Staleness
❌ **Wrong**: Forget to call `update_grid_pressures()` between steps
✅ **Right**: Call after mass-changing operations (flux)

Stale pressure → incorrect flux calculations.

### 4. CFL Violation
**Symptom**: Simulation explodes (NaN, Inf)
**Check**: Log `dt` and `u_max` - should see many small substeps if velocities are high
**Debug**: `println!("dt={}, u_max={}", dt, u_max);`

### 5. Boundary Confusion
- **Walls**: No-slip (u=0, v=0), zero-gradient for scalars
- **Open**: Vacuum (ρ=0) or zero-gradient (closed room)
- Make sure `collision_map` matches intended physics

---

## Performance

**Typical frame** (60 FPS, 31×16 grid):
- Substeps: 3-10 (velocity-dependent)
- Per substep:
  - Advection: ~0.5ms (bilinear interpolation)
  - Diffusion: ~0.2ms (5-point stencil)
  - Compression: ~0.1ms (divergence sampling)
- **Total: 2-5ms per frame** (within 16.7ms budget)

**Bottleneck**: Advection (bilinear interpolation pressure on cache)
**Optimization ideas**: cache-friendly tiling, SIMD interpolation, reduce substeps where possible

**Memory**: ~110 KB total (7 grid copies × 500 cells × 32 bytes/cell)

---

## Testing

Run step tests:
```bash
cargo test atmosphere::steps
```

**Existing tests**:
- `advection.rs`: MacCormack vs semi-Lagrangian, bilinear interpolation
- `compression.rs`: Expansion cooling and compression heating sanity checks
- `grid.rs`: Pressure calculation sanity check

**Recommended additional tests**:
- Vacuum-to-atmosphere interface
- Boundary condition correctness
- Conservation over multiple timesteps
- Divergence monitoring thresholds

---

## Constants Reference

From `../constants.rs`:

```rust
// Diffusion
KINEMATIC_VISCOSITY: 1.5e-5      // m²/s
THERMAL_DIFFUSIVITY: 2.1e-5      // m²/s
GAS_DIFFUSIVITY: 5.0e-4          // m²/s

// Thermodynamics
GAMMA: 1.4                      // Cp/Cv ratio

// CFL
CFL_NUMBER: 0.5
```

Legacy constants for the removed pressure flux (`PRESSURE_FLUX_CONDUCTANCE`, `PRESSURE_FLUX_MAX_FRACTION`) and Poisson projection (`POISSON_MAX_ITERATIONS`, `POISSON_TOLERANCE`, `SOR_OMEGA`) have been deleted from the codebase; refer to commit history if you need their last-known values.

---

## Further Reading

**Fluid Simulation**:
- Bridson, *Fluid Simulation for Computer Graphics* (2015)
- Stam, *Stable Fluids* (SIGGRAPH 1999)
- Fedkiw, *Visual Simulation of Smoke* (SIGGRAPH 2001)

**Numerical Methods**:
- LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002)
- Chorin, *Numerical Solution of the Navier-Stokes Equations* (1968)

**Operator Splitting**:
- Strang splitting (symmetric, 2nd-order)
- Godunov splitting (our approach, 1st-order, simpler)

---

## Related Files

- `../grid.rs` - `AtmosphereGrid` and `AtmosphereCell` definitions
- `../simulation.rs` - Main loop calling these steps
- `../constants.rs` - Physical constants
- `../monitoring.rs` - Conservation tracking
- `../../tilemap.rs` - `TileCollisionMap` for walls
