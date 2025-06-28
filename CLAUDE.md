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

## Important Notes

- No build system (build.zig) exists yet - will need to be created when implementing Zig code
- The `flake.nix` placeholder package (`pkgs.hello`) should be replaced with actual Zig package
- Development environment loads from `../nvim#zig`, suggesting a shared Neovim+Zig setup