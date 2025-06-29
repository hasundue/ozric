# Ozric

An Ornstein-Zernike (OZ) equation solver using ceres-solver, cross-compiled with Zig.

## Overview

Ozric leverages the ceres-solver optimization library to solve nonlinear systems arising from Ornstein-Zernike equations in statistical mechanics and liquid theory. The project uses Zig as the integration layer for cross-platform compilation and C++ library binding.

## Quick Start

### Using Nix (Recommended)

```bash
# Enter development environment
nix develop

# Build and run
nix build
nix run .

# Run tests
nix run .#test
```

### Using Zig directly

```bash
# Build the project
zig build

# Run the executable
zig build run

# Run tests
zig build test
```

## Architecture

- **Core**: Zig application with dual library/executable structure
- **Optimization**: ceres-solver integration for nonlinear system solving
- **Build system**: zig2nix for cross-platform Nix integration
- **Development**: Automated formatting, pre-commit hooks, comprehensive testing

## Status

Foundation implemented - ready for ceres-solver integration and OZ equation system development.
