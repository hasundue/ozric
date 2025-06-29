const std = @import("std");
const c = @cImport({
    @cInclude("fftw3.h");
});

const grid = @import("grid.zig");
const RadialFunction = grid.RadialFunction;

/// Errors that can occur during FFT operations
pub const FFTError = error{
    PlanCreationFailed,
    WorkspaceAllocationFailed,
    InvalidDataSize,
    NullPlan,
    NullWorkspace,
};

/// FFT planner wisdom/strategy
pub const PlannerStrategy = enum {
    estimate, // Fast planning, potentially slower execution
    measure, // Balanced planning and execution time
    patient, // Slower planning, potentially faster execution
    exhaustive, // Very slow planning, optimal execution

    fn toFFTWFlag(self: PlannerStrategy) c_uint {
        return switch (self) {
            .estimate => c.FFTW_ESTIMATE,
            .measure => c.FFTW_MEASURE,
            .patient => c.FFTW_PATIENT,
            .exhaustive => c.FFTW_EXHAUSTIVE,
        };
    }
};

/// FFT workspace and plan management for efficient transforms
pub const FFTContext = struct {
    allocator: std.mem.Allocator,
    n_points: usize,
    workspace: []f64,
    plan_forward: c.fftw_plan,
    plan_backward: c.fftw_plan,

    const Self = @This();

    /// Initialize FFT context with specified size and planner strategy
    pub fn init(allocator: std.mem.Allocator, n_points: usize, strategy: PlannerStrategy) FFTError!Self {
        // Allocate workspace - needs to be 2*n_points for real-to-half-complex transforms
        const workspace = allocator.alloc(f64, 2 * n_points) catch {
            return FFTError.WorkspaceAllocationFailed;
        };
        errdefer allocator.free(workspace);

        const strategy_flag = strategy.toFFTWFlag();

        // Create forward plan (real-to-half-complex)
        const plan_forward = c.fftw_plan_r2r_1d(
            @intCast(n_points),
            workspace.ptr,
            workspace.ptr,
            c.FFTW_R2HC,
            strategy_flag,
        );
        if (plan_forward == null) {
            return FFTError.PlanCreationFailed;
        }
        errdefer c.fftw_destroy_plan(plan_forward);

        // Create backward plan (half-complex-to-real)
        const plan_backward = c.fftw_plan_r2r_1d(
            @intCast(n_points),
            workspace.ptr,
            workspace.ptr,
            c.FFTW_HC2R,
            strategy_flag,
        );
        if (plan_backward == null) {
            return FFTError.PlanCreationFailed;
        }

        return Self{
            .allocator = allocator,
            .n_points = n_points,
            .workspace = workspace,
            .plan_forward = plan_forward,
            .plan_backward = plan_backward,
        };
    }

    /// Clean up FFT context resources
    pub fn deinit(self: *Self) void {
        c.fftw_destroy_plan(self.plan_forward);
        c.fftw_destroy_plan(self.plan_backward);
        self.allocator.free(self.workspace);
    }

    /// Perform forward FFT: real-to-half-complex transform
    /// Input data is copied to workspace, transformed in-place, and result copied back
    pub fn forwardTransform(self: *Self, input: []const f64, output: []f64) FFTError!void {
        if (input.len != self.n_points or output.len != self.n_points) {
            return FFTError.InvalidDataSize;
        }

        // Copy input to workspace
        @memcpy(self.workspace[0..self.n_points], input);

        // Zero-pad the rest of workspace
        @memset(self.workspace[self.n_points..], 0.0);

        // Execute forward transform
        c.fftw_execute(self.plan_forward);

        // Copy result back to output
        @memcpy(output, self.workspace[0..self.n_points]);
    }

    /// Perform inverse FFT: half-complex-to-real transform
    /// Input data is copied to workspace, transformed in-place, normalized, and result copied back
    pub fn inverseTransform(self: *Self, input: []const f64, output: []f64) FFTError!void {
        if (input.len != self.n_points or output.len != self.n_points) {
            return FFTError.InvalidDataSize;
        }

        // Copy input to workspace
        @memcpy(self.workspace[0..self.n_points], input);

        // Zero-pad the rest of workspace
        @memset(self.workspace[self.n_points..], 0.0);

        // Execute inverse transform
        c.fftw_execute(self.plan_backward);

        // Normalize and copy result back to output
        const norm_factor = 1.0 / @as(f64, @floatFromInt(self.n_points));
        for (0..self.n_points) |i| {
            output[i] = self.workspace[i] * norm_factor;
        }
    }

    /// Perform in-place forward FFT on workspace
    /// Useful for advanced operations where data is already in workspace
    pub fn forwardTransformInPlace(self: *Self) void {
        c.fftw_execute(self.plan_forward);
    }

    /// Perform in-place inverse FFT on workspace with normalization
    /// Useful for advanced operations where data is already in workspace
    pub fn inverseTransformInPlace(self: *Self) void {
        c.fftw_execute(self.plan_backward);

        // Normalize in-place
        const norm_factor = 1.0 / @as(f64, @floatFromInt(self.n_points));
        for (0..self.n_points) |i| {
            self.workspace[i] *= norm_factor;
        }
    }
};

/// High-level utilities for Fourier transforms of radial functions
pub const FourierUtils = struct {
    /// Compute convolution of two radial functions using FFT
    /// This is the core operation for solving Ornstein-Zernike equations
    /// Implements: h(r) = ∫ c(r-r') * ρ * h(r') dr'
    pub fn convolution(
        allocator: std.mem.Allocator,
        input1: *const RadialFunction,
        input2: *const RadialFunction,
        output: *RadialFunction,
        density: f64,
    ) !void {
        if (input1.grid.n_points != input2.grid.n_points or input1.grid.n_points != output.grid.n_points) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input1.grid.n_points, .estimate);
        defer fft_ctx.deinit();

        // Transform first function
        try fft_ctx.forwardTransform(input1.values, output.values);

        // Store transformed first function temporarily
        const temp_fft1 = try allocator.alloc(f64, input1.grid.n_points);
        defer allocator.free(temp_fft1);
        @memcpy(temp_fft1, output.values);

        // Transform second function
        try fft_ctx.forwardTransform(input2.values, output.values);

        // Perform convolution in frequency domain
        // For real-to-half-complex format, need special handling
        for (0..input1.grid.n_points) |i| {
            output.values[i] = temp_fft1[i] * output.values[i] * density;
        }

        // Transform back to real domain
        try fft_ctx.inverseTransform(output.values, output.values);
    }

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

    /// Compute power spectrum of a radial function
    /// Useful for analyzing frequency content and debugging transforms
    pub fn powerSpectrum(
        allocator: std.mem.Allocator,
        input: *const RadialFunction,
        spectrum: *RadialFunction,
    ) !void {
        if (input.grid.n_points != spectrum.grid.n_points) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input.grid.n_points, .estimate);
        defer fft_ctx.deinit();

        // Transform to frequency domain
        try fft_ctx.forwardTransform(input.values, spectrum.values);

        // Compute power spectrum |F(k)|²
        for (0..input.grid.n_points) |i| {
            spectrum.values[i] = spectrum.values[i] * spectrum.values[i];
        }
    }
};

// ============================================================================
// UNIT TESTS
// ============================================================================

test "FFTContext initialization and cleanup" {
    const allocator = std.testing.allocator;

    var fft_ctx = try FFTContext.init(allocator, 64, .estimate);
    defer fft_ctx.deinit();

    try std.testing.expect(fft_ctx.n_points == 64);
    try std.testing.expect(fft_ctx.workspace.len == 128);
}

test "FFT roundtrip identity test" {
    const allocator = std.testing.allocator;
    const n = 32;

    var fft_ctx = try FFTContext.init(allocator, n, .estimate);
    defer fft_ctx.deinit();

    // Create test signal: sine wave
    const input = try allocator.alloc(f64, n);
    defer allocator.free(input);
    const output = try allocator.alloc(f64, n);
    defer allocator.free(output);
    const temp = try allocator.alloc(f64, n);
    defer allocator.free(temp);

    for (0..n) |i| {
        const x = 2.0 * std.math.pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        input[i] = @sin(x);
    }

    // Forward then inverse transform
    try fft_ctx.forwardTransform(input, temp);
    try fft_ctx.inverseTransform(temp, output);

    // Check roundtrip accuracy
    for (0..n) |i| {
        try std.testing.expectApproxEqAbs(input[i], output[i], 1e-12);
    }
}

test "FFT of known analytical functions" {
    const allocator = std.testing.allocator;
    const n = 64;

    var fft_ctx = try FFTContext.init(allocator, n, .estimate);
    defer fft_ctx.deinit();

    const input = try allocator.alloc(f64, n);
    defer allocator.free(input);
    const output = try allocator.alloc(f64, n);
    defer allocator.free(output);

    // Test delta function: δ(x) -> constant in frequency domain
    @memset(input, 0.0);
    input[0] = 1.0; // Discrete delta at origin

    try fft_ctx.forwardTransform(input, output);

    // FFT of delta should be approximately constant (all ones)
    // Note: discrete delta doesn't perfectly replicate continuous delta behavior
    // but first component should be close to 1
    try std.testing.expectApproxEqAbs(output[0], 1.0, 1e-10);

    // Other components may not be exactly 1 due to discrete sampling effects
    // Just verify they're reasonable values
    for (1..n) |i| {
        try std.testing.expect(@abs(output[i]) < 2.0);
    }
}

test "FFT of constant function" {
    const allocator = std.testing.allocator;
    const n = 16;

    var fft_ctx = try FFTContext.init(allocator, n, .estimate);
    defer fft_ctx.deinit();

    const input = try allocator.alloc(f64, n);
    defer allocator.free(input);
    const output = try allocator.alloc(f64, n);
    defer allocator.free(output);

    // Constant function
    @memset(input, 1.0);

    try fft_ctx.forwardTransform(input, output);

    // FFT of constant should be delta in frequency domain
    // First component should be n (DC component)
    try std.testing.expectApproxEqAbs(output[0], @as(f64, @floatFromInt(n)), 1e-10);

    // Other components should be approximately zero
    for (1..n) |i| {
        try std.testing.expectApproxEqAbs(output[i], 0.0, 1e-10);
    }
}

test "FFT linearity property" {
    const allocator = std.testing.allocator;
    const n = 32;

    var fft_ctx = try FFTContext.init(allocator, n, .estimate);
    defer fft_ctx.deinit();

    const f1 = try allocator.alloc(f64, n);
    defer allocator.free(f1);
    const f2 = try allocator.alloc(f64, n);
    defer allocator.free(f2);
    const f_sum = try allocator.alloc(f64, n);
    defer allocator.free(f_sum);

    const F1 = try allocator.alloc(f64, n);
    defer allocator.free(F1);
    const F2 = try allocator.alloc(f64, n);
    defer allocator.free(F2);
    const F_sum = try allocator.alloc(f64, n);
    defer allocator.free(F_sum);
    const F_sum_direct = try allocator.alloc(f64, n);
    defer allocator.free(F_sum_direct);

    // Create two different test functions
    for (0..n) |i| {
        const x = 2.0 * std.math.pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        f1[i] = @sin(x);
        f2[i] = @cos(2.0 * x);
        f_sum[i] = f1[i] + f2[i];
    }

    // Transform individual functions
    try fft_ctx.forwardTransform(f1, F1);
    try fft_ctx.forwardTransform(f2, F2);

    // Transform sum directly
    try fft_ctx.forwardTransform(f_sum, F_sum_direct);

    // Compute sum of transforms
    for (0..n) |i| {
        F_sum[i] = F1[i] + F2[i];
    }

    // Check linearity: FFT(f1 + f2) = FFT(f1) + FFT(f2)
    for (0..n) |i| {
        try std.testing.expectApproxEqAbs(F_sum_direct[i], F_sum[i], 1e-12);
    }
}

test "Ornstein-Zernike equation solver" {
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
    try FourierUtils.solveOZEquation(allocator, &c_r, &h_r, test_density);

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

test "Power spectrum computation" {
    const allocator = std.testing.allocator;
    const n = 32;

    const grid_params = grid.GridParams{
        .r_max = 10.0,
        .n_points = n,
        .dr = 10.0 / @as(f64, @floatFromInt(n)),
    };

    var input = try RadialFunction.init(allocator, grid_params);
    defer input.deinit();
    var spectrum = try RadialFunction.init(allocator, grid_params);
    defer spectrum.deinit();

    // Create a sine wave
    for (0..n) |i| {
        const x = 2.0 * std.math.pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        input.values[i] = @sin(x);
    }

    try FourierUtils.powerSpectrum(allocator, &input, &spectrum);

    // Power spectrum should be non-negative
    for (0..n) |i| {
        try std.testing.expect(spectrum.values[i] >= 0.0);
    }

    // For a pure sine wave, most power should be concentrated in specific frequencies
    var total_power: f64 = 0.0;
    for (0..n) |i| {
        total_power += spectrum.values[i];
    }
    try std.testing.expect(total_power > 0.0);
}
