# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Ozric is an Ornstein-Zernike (OZ) equation solver that leverages ceres-solver for optimization, cross-compiled with Zig. The project **statically compiles ceres-solver from source** using Zig's build system, integrating C++ numerical optimization libraries for solving nonlinear systems in statistical mechanics and liquid theory.

### Core Design Principle
**IMPORTANT**: This project compiles ceres-solver from source using Zig, NOT using pre-built libraries. Static compilation with Zig is fundamental to the project architecture.

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

### Ceres-Solver Integration Status ✅

The project has successfully achieved **complete ceres-solver 2.2.0 integration** with dual-target support:

#### **Native Build (Full-Featured)**
- ✅ **Complete Ceres 2.2.0 integration** with all template specializations
- ✅ **Threading support** with full optimization capabilities  
- ✅ **159MB executable** with comprehensive functionality
- ✅ **Working hello world** demonstrating successful optimization

#### **WASM Build (Browser-Compatible)**
- ✅ **4.0MB WASM binary** for web deployment
- ✅ **Ultra-minimal Ceres subset** avoiding threading dependencies
- ✅ **Basic functionality** (version detection, cost/loss functions)
- ✅ **Browser-ready** single artifact deployment

#### **Threading Architecture Limitations**

**CRITICAL**: Ceres 2.2.0 has **mandatory threading dependencies** that cannot be disabled:

```cpp
// ceres/problem_impl.h:49 - EVERY Problem includes threading
#include "ceres/context_impl.h"

// ceres/context_impl.h:68 - Threading always present  
ThreadPool thread_pool;  // Unconditional member

// ceres/thread_pool.h:35-36 - Headers require threading primitives
#include <mutex>
#include <thread>
```

**Why `CERES_NO_THREADS` doesn't work:**
- ❌ **Completely removed in Ceres 2.2.0** - flag is a no-op
- ❌ **Zero conditional compilation** - no `#ifdef CERES_NO_THREADS` blocks exist
- ❌ **Architectural change** - threading is now mandatory, not optional
- ❌ **Documentation lag** - old tutorials still reference the defunct flag

From Ceres 2.2.0 changelog:
> "OpenMP and NO_THREADING backends have been removed. C++ threads is how all threading is done."

This means **full Ceres functionality in WASM requires threading polyfills** or accepting larger binary sizes.

### Build Architecture

#### **Dual-Target Build System**
- **Native**: Full ceres source set (~80 files) with template specializations
- **WASM**: Minimal source set (6 files) avoiding threading dependencies  
- **Shared**: Common C++ wrapper (`src/solver.cc`, `include/ceres.h`)
- **Module separation**: Prevents threading contamination between targets

#### **Template Specialization Management**
- **Programmatic generation** of 38 template specialization files
- **Native-only**: SchurEliminator and PartitionedMatrixView specializations
- **WASM exclusion**: Avoids template-induced threading dependencies

### Integration Achievements
- ✅ **Static compilation** of ceres-solver from source with Zig
- ✅ **Cross-platform support** (native + WASM)
- ✅ **Zero external dependencies** - fully self-contained builds
- ✅ **Working optimization** - hello world converges correctly
- ✅ **Browser deployment ready** - single 4MB WASM file