# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive observability infrastructure for atmosphere simulation
  - Real-time metrics collection and monitoring system
  - CSV export functionality for simulation data
  - Visualization tools for pressure and velocity fields
  - Detailed diagnostics with mass conservation and divergence tracking
- Integration tests for atmosphere simulation (`tests/atmosphere_integration.rs`)
- Fluid dynamics fundamental verification tests (`tests/fluid_dynamics_fundamentals.rs`)
- Detailed documentation for atmosphere simulation fixes and bug analysis
- CONTRIBUTING.md with development workflow and commit conventions
- AGENTS.md for AI agent collaboration guidelines

### Changed
- Refactored advection step with improved numerical stability
- Enhanced compression step implementation
- Updated atmosphere simulation constants for better physical accuracy
- Improved monitoring and debug output formatting
- Updated dependencies in Cargo.toml

### Removed
- Legacy pressure projection implementation (`pressure_projection.rs`)
- Legacy pressure flux implementation (`pressure_flux.rs`)

### Fixed
- Mass conservation issues in atmosphere simulation
- Momentum transfer calculations in advection step
