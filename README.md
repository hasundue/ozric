# Ozric

A proof-of-concept solver for Ornstein-Zernike equations implemented in Zig.

## Overview

Ozric solves the Ornstein-Zernike equation for liquid theory and statistical mechanics:

```
h(r) = c(r) + ρ ∫ c(|r-r'|) h(r') dr'
```

Where:
- `h(r)` is the total correlation function
- `c(r)` is the direct correlation function  
- `ρ` is the number density
- The integral represents convolution over all space

## Features

### Closure Relations
- **Hypernetted Chain (HNC)**: `c(r) = h(r) - ln(g(r)) - βu(r)`
- **Percus-Yevick (PY)**: `c(r) = (1 - exp(βu(r))) * g(r)`

### Interaction Potentials
- **Hard Sphere**: Step function potential for excluded volume
- **Lennard-Jones**: `u(r) = 4ε[(σ/r)¹² - (σ/r)⁶]` for realistic molecular interactions

### Solver Methods
- **Simple Convolution**: Educational implementation with direct integration
- **FFT-based**: Efficient Fourier-space solver using FFTW for convolution

### Convergence Acceleration
- **Simple Mixing**: Basic linear mixing `h_new = αh_new + (1-α)h_old`
- **Adaptive Mixing**: Dynamic mixing factor based on convergence behavior
- **Anderson Acceleration**: History-based extrapolation for faster convergence

### Numerical Stability
- NaN detection and recovery in iteration process
- Potential energy clamping to prevent overflow
- Conservative initial guesses with mean field theory
- Percus-Yevick analytical approximations for hard spheres

## Dependencies

- **Zig**: Modern systems programming language
- **BLAS**: Basic Linear Algebra Subprograms
- **LAPACK**: Linear Algebra Package  
- **FFTW**: Fast Fourier Transform library

## Usage

```bash
# Build the project
zig build

# Run tests
zig build test

# Run the executable
zig build run
```

## Development

Uses Nix flakes for reproducible development environment:

```bash
# Enter development shell
nix develop

# Format code
treefmt
```

## Project Structure

```
src/
├── main.zig          # Executable entry point
├── root.zig          # Public API re-exports and tests
├── potentials.zig    # Interaction potential interface and implementations
├── grid.zig          # Grid parameters and radial function data structures
├── convergence.zig   # Convergence acceleration algorithms
└── solver.zig        # Main OZ equation solver with closure relations
tests/
├── integration_test.zig  # Integration tests
build.zig             # Build configuration with C library linking
```

### Module Organization

- **`potentials.zig`**: Extensible potential interface with Hard Sphere and Lennard-Jones implementations
- **`grid.zig`**: Radial grid discretization, correlation functions, and BLAS-integrated vectors  
- **`convergence.zig`**: Convergence acceleration methods (simple mixing, adaptive mixing, Anderson)
- **`solver.zig`**: Core OZ equation solver with closure relations and FFT-based algorithms
- **`root.zig`**: Clean public API that re-exports all modules for easy consumption

## Mathematical Background

The implementation solves the Ornstein-Zernike equation coupled with closure relations to determine the radial distribution function `g(r) = h(r) + 1`, which describes the probability of finding particles at distance `r` in liquids and dense gases.

This is fundamental for:
- Liquid state theory
- Statistical mechanics
- Molecular dynamics validation
- Equation of state calculations
