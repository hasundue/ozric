# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ozric is an Ornstein-Zernike (OZ) equation solver that leverages ceres-solver for optimization, cross-compiled with Zig. The project integrates C++ numerical optimization libraries for solving nonlinear systems in statistical mechanics and liquid theory.

## Development Environment

### Nix Flakes Setup
This project uses Nix flakes for reproducible development environments:

```bash
# Enter development shell (includes ZLS, treefmt, nil)
nix develop

# Build and run the project
nix build
nix run .

# Format code
treefmt
```

### direnv Integration
The repository uses direnv with two flakes:
- `../nvim#zig` - Neovim with Zig language support
- `.` - Local project dependencies

### Code Formatting
- Pre-commit hooks automatically run `treefmt` on commit
- Currently configured for Nix files; Zig formatting can be added to treefmt.nix
- Use `treefmt --fail-on-change --no-cache` to check formatting
- Git hooks are managed by git-hooks.nix (do not manually edit .pre-commit-config.yaml)

## Build System

### Zig Build Configuration
The project uses Zig's native build system with dual-module architecture:

```bash
# Build the project
zig build

# Run the executable
zig build run

# Run all tests (library + executable)
zig build test

# Generate documentation
nix run .#docs
```

### Build Targets and Structure
- **Library module**: `src/root.zig` - Core functionality with exported functions
- **Executable module**: `src/main.zig` - Application entry point
- **Module dependency**: Executable imports library as "default_lib"
- **Test coverage**: Both modules include comprehensive unit tests

### Nix Integration
```bash
# Build with Nix (two variants available)
nix build          # Default Nix-friendly binaries
nix build .#foreign # Musl-linked binaries for non-Nix distribution

# Development and testing
nix run .#build    # Build with arguments
nix run .#test     # Run tests via Nix
```

## Architecture

### Current Implementation
- **Language**: Zig (minimum version 0.14.1)
- **Build system**: Complete zig2nix integration with cross-platform support
- **Code structure**: Library/executable hybrid with clean module separation
- **Testing infrastructure**: Unit tests, integration tests, and fuzz testing examples
- **Package management**: Ready for external dependencies via build.zig.zon

### Code Quality Pipeline
- **Automated formatting**: treefmt with pre-commit hooks
- **Memory safety**: Tests include memory leak detection
- **Cross-compilation**: Supports all Zig target platforms via Nix
- **IDE support**: ZLS (Zig Language Server) and Neovim integration

### Future Integration Points
- **ceres-solver integration**: Architecture prepared for C++ optimization library binding
- **Cross-compilation**: Zig build system ready for C++ library linking
- **Domain-specific code**: Foundation prepared for OZ equation system implementation
- **Package naming**: Currently uses template names ("default") - ready for project-specific naming