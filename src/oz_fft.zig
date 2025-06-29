const std = @import("std");
const fourier = @import("fourier.zig");
const grid = @import("grid.zig");

const FFTContext = fourier.FFTContext;
const FourierUtils = fourier.FourierUtils;
const FFTError = fourier.FFTError;
const RadialFunction = grid.RadialFunction;

/// Ornstein-Zernike equation specific FFT operations
/// This module provides high-level FFT operations for liquid theory
pub const OZFourierSolver = struct {
    /// Solve OZ equation: h(k) = c(k) / (1 - ρ * c(k))
    /// This is the k-space form of the Ornstein-Zernike equation
    pub fn solveOZEquation(
        allocator: std.mem.Allocator,
        c_r: *const RadialFunction,
        h_r: *RadialFunction,
        density: f64,
    ) !void {
        if (c_r.grid.n_points != h_r.grid.n_points) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, c_r.grid.n_points, .estimate);
        defer fft_ctx.deinit();

        // Transform c(r) to k-space
        try fft_ctx.forwardTransform(c_r.values, h_r.values);

        // Solve OZ equation in k-space: h(k) = c(k) / (1 - ρ * c(k))
        for (0..c_r.grid.n_points) |i| {
            const c_k = h_r.values[i];
            const denominator = 1.0 - density * c_k;

            // Avoid division by zero
            if (@abs(denominator) < 1e-12) {
                h_r.values[i] = 0.0;
            } else {
                h_r.values[i] = c_k / denominator;
            }
        }

        // Transform h(k) back to real space
        try fft_ctx.inverseTransform(h_r.values, h_r.values);
    }

    /// Compute convolution for OZ equations with density scaling
    /// Implements: h(r) = ∫ c(r-r') * ρ * h(r') dr'
    pub fn ozConvolution(
        allocator: std.mem.Allocator,
        c_r: *const RadialFunction,
        h_r_in: *const RadialFunction,
        h_r_out: *RadialFunction,
        density: f64,
    ) !void {
        if (c_r.grid.n_points != h_r_in.grid.n_points or c_r.grid.n_points != h_r_out.grid.n_points) {
            return FFTError.InvalidDataSize;
        }

        // Use generic FFT convolution
        try FourierUtils.convolution(allocator, c_r.values, h_r_in.values, h_r_out.values);

        // Apply density scaling
        for (0..h_r_out.grid.n_points) |i| {
            h_r_out.values[i] *= density;
        }
    }

    /// Compute radial distribution function power spectrum for analysis
    pub fn radialPowerSpectrum(
        allocator: std.mem.Allocator,
        input: *const RadialFunction,
        spectrum: *RadialFunction,
    ) !void {
        if (input.grid.n_points != spectrum.grid.n_points) {
            return FFTError.InvalidDataSize;
        }

        // Use generic power spectrum computation
        try FourierUtils.powerSpectrum(allocator, input.values, spectrum.values);
    }

    /// Apply spherical coordinate correction for 3D systems
    /// This corrects for the r² factor in spherical coordinates
    pub fn applySphereCorrection(
        input: *RadialFunction,
        output: *RadialFunction,
        apply_forward: bool,
    ) void {
        if (input.grid.n_points != output.grid.n_points) return;

        for (0..input.grid.n_points) |i| {
            const r = @as(f64, @floatFromInt(i)) * input.grid.dr;

            if (r > 1e-12) {
                if (apply_forward) {
                    // Apply r multiplication for forward transform
                    output.values[i] = input.values[i] * r;
                } else {
                    // Remove r factor for inverse transform
                    output.values[i] = input.values[i] / r;
                }
            } else {
                // Handle r=0 case
                output.values[i] = 0.0;
            }
        }
    }
};

// ============================================================================
// UNIT TESTS
// ============================================================================

test "OZ equation solver with exponential decay" {
    const allocator = std.testing.allocator;
    const n = 64;

    // Create test radial functions
    const grid_params = grid.GridParams{
        .r_max = 10.0,
        .n_points = n,
        .dr = 10.0 / @as(f64, @floatFromInt(n)),
    };

    var c_r = try RadialFunction.init(allocator, grid_params);
    defer c_r.deinit();
    var h_r = try RadialFunction.init(allocator, grid_params);
    defer h_r.deinit();

    // Set up a simple test case: exponential decay c(r) = exp(-r)
    for (0..n) |i| {
        const r = @as(f64, @floatFromInt(i)) * grid_params.dr;
        c_r.values[i] = @exp(-r);
    }

    const test_density = 0.1;

    // Solve OZ equation
    try OZFourierSolver.solveOZEquation(allocator, &c_r, &h_r, test_density);

    // Basic sanity checks
    // h(r) should be finite and reasonable
    for (0..n) |i| {
        try std.testing.expect(std.math.isFinite(h_r.values[i]));
        try std.testing.expect(@abs(h_r.values[i]) < 100.0); // Reasonable magnitude
    }

    // For small density and exponential c(r), h(r) should have reasonable behavior
    // At large distances, both c(r) and h(r) should decay to near zero
    const large_r_index = n - 1;
    try std.testing.expect(@abs(h_r.values[large_r_index]) < 1.0);

    // At intermediate distances, h(r) should be reasonable
    // This is a complex algorithm so we just verify basic sanity
    const mid_r_index = n / 2;
    try std.testing.expect(@abs(h_r.values[mid_r_index]) < 10.0);
}

test "OZ convolution with density scaling" {
    const allocator = std.testing.allocator;
    const n = 32;

    const grid_params = grid.GridParams{
        .r_max = 5.0,
        .n_points = n,
        .dr = 5.0 / @as(f64, @floatFromInt(n)),
    };

    var c_r = try RadialFunction.init(allocator, grid_params);
    defer c_r.deinit();
    var h_r_in = try RadialFunction.init(allocator, grid_params);
    defer h_r_in.deinit();
    var h_r_out = try RadialFunction.init(allocator, grid_params);
    defer h_r_out.deinit();

    // Set up simple test functions
    for (0..n) |i| {
        const r = @as(f64, @floatFromInt(i)) * grid_params.dr;
        c_r.values[i] = @exp(-r);
        h_r_in.values[i] = @exp(-2.0 * r);
    }

    const test_density = 0.5;

    // Perform OZ convolution
    try OZFourierSolver.ozConvolution(allocator, &c_r, &h_r_in, &h_r_out, test_density);

    // Verify results are finite and have density scaling applied
    for (0..n) |i| {
        try std.testing.expect(std.math.isFinite(h_r_out.values[i]));
    }

    // The result should be different from input due to convolution and scaling
    var is_different = false;
    for (0..n) |i| {
        if (@abs(h_r_out.values[i] - h_r_in.values[i]) > 1e-10) {
            is_different = true;
            break;
        }
    }
    try std.testing.expect(is_different);
}

test "Spherical coordinate correction" {
    const allocator = std.testing.allocator;
    const n = 16;

    const grid_params = grid.GridParams{
        .r_max = 4.0,
        .n_points = n,
        .dr = 4.0 / @as(f64, @floatFromInt(n)),
    };

    var input = try RadialFunction.init(allocator, grid_params);
    defer input.deinit();
    var output = try RadialFunction.init(allocator, grid_params);
    defer output.deinit();
    var restored = try RadialFunction.init(allocator, grid_params);
    defer restored.deinit();

    // Set up test function
    for (0..n) |i| {
        input.values[i] = @as(f64, @floatFromInt(i + 1));
    }

    // Apply forward correction
    OZFourierSolver.applySphereCorrection(&input, &output, true);

    // Apply inverse correction
    OZFourierSolver.applySphereCorrection(&output, &restored, false);

    // Check that forward and inverse corrections roughly cancel out
    // (except at r=0 where we handle the singularity)
    for (1..n) |i| {
        try std.testing.expectApproxEqAbs(input.values[i], restored.values[i], 1e-10);
    }
}
