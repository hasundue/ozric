# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ozric is a proof-of-concept solver for Ornstein-Zernike equations implemented in Zig. The project integrates C mathematics libraries (BLAS, LAPACK, FFTW, nlopt) for numerical computations in statistical mechanics and liquid theory. All C libraries are statically linked for WebAssembly compatibility.

## Development Environment

### Nix Flakes Setup
This project uses Nix flakes for reproducible development environments:

```bash
# Enter development shell
nix develop

# Build the project - MUST use nix develop --command for proper library linking
nix develop --command zig build

# Run the executable
nix develop --command zig build run

# Run tests
nix develop --command zig build test

# Format code
treefmt
```

**IMPORTANT**: Due to static linking of C libraries (BLAS, LAPACK, FFTW, nlopt), all Zig build commands must be run with `nix develop --command` to ensure proper pkg-config environment setup.

### direnv Integration
The repository uses direnv with two flakes:
- `../nvim#zig` - Neovim with Zig language support
- `.` - Local project dependencies

### Code Formatting
- Pre-commit hooks automatically run `treefmt` on commit
- Currently configured for Nix files; will need Zig formatting when source code is added
- Use `treefmt --fail-on-change --no-cache` to check formatting

## Architecture

### Current State
The project includes:
- Cross-platform Nix development environment (Linux x86_64/aarch64, macOS x86_64/aarch64)  
- Automated code formatting and git hooks
- MIT license (Copyright 2025 Shun Ueda)
- Working Zig implementation with C library integration
- Static linking setup for WebAssembly compatibility

### Architecture
- **Language**: Zig for the main solver implementation
- **Dependencies**: C mathematics libraries (BLAS, LAPACK, FFTW, nlopt) - all statically linked
- **Domain**: Ornstein-Zernike equation solving for statistical mechanics
- **Features**: Vector operations (BLAS), linear algebra (LAPACK), Fourier transforms (FFTW), nonlinear optimization (nlopt)

## Build Commands

**IMPORTANT**: All Zig build commands must be prefixed with `nix develop --command` for proper library linking.

```bash
# Build the project
nix develop --command zig build

# Run the executable
nix develop --command zig build run

# Run tests
nix develop --command zig build test

# Build in release mode
nix develop --command zig build -Doptimize=ReleaseFast

# Clean build artifacts
rm -rf zig-cache zig-out
```

## Project Structure

```
src/
├── main.zig          # Executable entry point
├── root.zig          # Library root with core functionality
tests/
├── integration_test.zig  # Integration tests
build.zig             # Build configuration
```

## Important Notes

- Project uses C mathematics libraries: BLAS, LAPACK, FFTW3, nlopt (all statically linked)
- Static linking enables WebAssembly compilation for the solver
- Tests include both unit tests (in src files) and integration tests (tests/ directory)  
- Development environment loads from `../nvim#zig`, suggesting a shared Neovim+Zig setup
- **CRITICAL**: Always use `nix develop --command` prefix for Zig builds to ensure proper pkg-config environment