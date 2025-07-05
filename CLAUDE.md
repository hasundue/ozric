# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ozric is a classical density functional theory (DFT) solver that leverages ceres-solver for optimization, cross-compiled with Zig. The project **statically compiles ceres-solver from source** using Zig's build system, integrating C++ numerical optimization libraries for solving variational problems in statistical mechanics, particularly free energy functionals and density distributions.

### Core Design Principle
**IMPORTANT**: This project compiles ceres-solver from source using Zig, NOT using pre-built libraries. Static compilation with Zig is fundamental to the project architecture.

## Development Environment

### Nix Flakes Setup
This project uses Nix flakes for reproducible development environments:

```bash
# Enter development shell (includes ZLS, clangd, treefmt)
nix develop

# Build the library
nix build
nix build .#foreign  # Musl-linked binaries for non-Nix distribution

# Format code
treefmt
```

### direnv Integration
The repository uses direnv with two flakes:
- `../nvim#zig` - Neovim with Zig language support
- `.` - Local project dependencies

### Code Formatting
- **Automated formatting**: treefmt with pre-commit hooks
- **Languages supported**: Zig (zig fmt), C++ (clang-format), Nix (nixfmt)
- **C++ style**: Google style with 4-space indentation (configured in `.clang-format`)
- **Zig style**: Standard zig fmt with 4-space indentation (non-configurable)
- **Consistency**: Both Zig and C++ use 4-space indentation for uniformity
- **Usage**: Run `treefmt` manually or automatic via pre-commit hooks
- **Validation**: Use `treefmt --fail-on-change --no-cache` to check formatting
- **Configuration**: Git hooks managed by git-hooks.nix (do not manually edit .pre-commit-config.yaml)

## Build System

### Zig Build Configuration
The project uses Zig's native build system as a **library-only** project:

```bash
# Build the library
zig build

# Run all tests (includes former main() functionality as test)
zig build test

# Generate compile_commands.json for clangd
zig build compile-commands

# Generate fish shell completions
zig build fish-completions

# Build WASM (use Nix for proper environment)
nix run .#build-wasm
```

### Build Targets and Structure
- **Library module**: `src/root.zig` - Core functionality with exported functions
- **C++ Integration**: `src/solver.cc` - C++ wrapper for ceres-solver functions
- **Build Module**: `src/build/ceres.zig` - Complex ceres-solver integration and WASM support
- **Headers**: `include/` - C interface headers for external usage
- **Test coverage**: Library includes comprehensive unit tests and former main() functionality

### Nix Integration
```bash
# Build with Nix (two variants available)
nix build          # Default Nix-friendly binaries
nix build .#foreign # Musl-linked binaries for non-Nix distribution

# Development and testing
nix run .#test      # Run tests via Nix
nix run .#build-wasm # Build WASM with proper Emscripten environment
```

## Architecture

### Current Implementation
- **Language**: Zig (minimum version 0.14.1)
- **Build system**: Complete zig2nix integration with cross-platform support
- **Code structure**: Library-only with C++ integration via ceres-solver
- **Testing infrastructure**: Unit tests, integration tests, and fuzz testing examples
- **Package management**: Ready for external dependencies via build.zig.zon

### Ceres-Solver Integration
- **Static compilation**: Builds ceres-solver 2.2.0 from source with minimal dependencies
- **Configuration**: Uses MINIGLOG, disables LAPACK, CUDA, SuiteSparse for portability
- **WASM support**: Full Emscripten compilation pipeline with sophisticated object file management
- **IDE support**: Automatic compile_commands.json generation for clangd C++ support

### Code Quality Pipeline
- **Automated formatting**: treefmt with pre-commit hooks
- **Memory safety**: Tests include memory leak detection
- **Cross-compilation**: Supports all Zig target platforms via Nix
- **IDE support**: ZLS (Zig Language Server) and clangd integration
