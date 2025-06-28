# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ozric is a proof-of-concept solver for Ornstein-Zernike equations implemented in Zig. The project integrates C mathematics libraries (BLAS, LAPACK, FFTW) for numerical computations in statistical mechanics and liquid theory.

## Development Environment

### Nix Flakes Setup
This project uses Nix flakes for reproducible development environments:

```bash
# Enter development shell
nix develop

# Build the project (once implemented)
nix build

# Format code
treefmt
```

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
The project is in initial setup phase with no Zig source code yet. The foundation includes:
- Cross-platform Nix development environment (Linux x86_64/aarch64, macOS x86_64/aarch64)
- Automated code formatting and git hooks
- MIT license (Copyright 2025 Shun Ueda)

### Planned Architecture
- **Language**: Zig for the main solver implementation
- **Dependencies**: C mathematics libraries (BLAS, LAPACK, FFTW)
- **Domain**: Ornstein-Zernike equation solving for statistical mechanics

## Build Commands

```bash
# Build the project
zig build

# Run the executable
zig build run

# Run tests
zig build test

# Build in release mode
zig build -Doptimize=ReleaseFast

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

- Project uses C mathematics libraries: BLAS, LAPACK, FFTW3
- Tests include both unit tests (in src files) and integration tests (tests/ directory)
- Development environment loads from `../nvim#zig`, suggesting a shared Neovim+Zig setup