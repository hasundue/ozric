# Ozric

A classical density functional theory (DFT) solver using ceres-solver, cross-compiled with Zig.

## Overview

Ozric leverages the ceres-solver optimization library to solve nonlinear systems arising from classical density functional theory. The project focuses on variational optimization problems in statistical mechanics, particularly free energy functionals and density distributions. Zig serves as the integration layer for cross-platform compilation and C++ library binding.

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
- **Optimization**: ceres-solver integration for variational DFT problems
- **Applications**: Free energy minimization, density profile optimization, phase transitions
- **Build system**: zig2nix for cross-platform Nix integration
- **Development**: Automated formatting, pre-commit hooks, comprehensive testing

## Development Tools

### IDE Support

**C++ Language Server (clangd)**:
- Automatically generates `compile_commands.json` when entering the development environment
- Provides full IDE features for C++ code including autocomplete, error checking, and go-to-definition

**Code Formatting**:
- Automated formatting with `treefmt` (Zig, C++, Nix)
- C++ formatting based on Google style with 4-space indentation (matches Zig)
- Pre-commit hooks automatically format code on commit
- Run manually: `treefmt`

**Fish Shell Completions**:
- Automatically generates shell completions for custom `zig build` steps
- To use: `source completions/zig.fish` in fish shell
- Available completions: `compile-commands`, `wasm`, `fish-completions`

### Available Build Steps

```bash
zig build                    # Default build
zig build run               # Run the application
zig build test              # Run tests
zig build wasm              # Build WebAssembly version
zig build compile-commands  # Generate clangd compilation database
zig build fish-completions  # Generate fish shell completions
treefmt                     # Format all code (Zig, C++, Nix)
```

## Status

Foundation implemented - ready for classical DFT functional development and variational optimization problems.
