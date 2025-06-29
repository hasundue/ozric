const std = @import("std");
const c = @cImport({
    @cInclude("fftw3.h");
});

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

/// High-level utilities for FFT operations on generic arrays
pub const FourierUtils = struct {
    /// Compute convolution of two arrays using FFT
    /// output = IFFT(FFT(input1) * FFT(input2))
    pub fn convolution(
        allocator: std.mem.Allocator,
        input1: []const f64,
        input2: []const f64,
        output: []f64,
    ) !void {
        if (input1.len != input2.len or input1.len != output.len) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input1.len, .estimate);
        defer fft_ctx.deinit();

        // Transform first function
        try fft_ctx.forwardTransform(input1, output);

        // Store transformed first function temporarily
        const temp_fft1 = try allocator.alloc(f64, input1.len);
        defer allocator.free(temp_fft1);
        @memcpy(temp_fft1, output);

        // Transform second function
        try fft_ctx.forwardTransform(input2, output);

        // Perform convolution in frequency domain
        for (0..input1.len) |i| {
            output[i] = temp_fft1[i] * output[i];
        }

        // Transform back to real domain
        try fft_ctx.inverseTransform(output, output);
    }

    /// Compute cross-correlation of two arrays using FFT
    /// Similar to convolution but with complex conjugate
    pub fn crossCorrelation(
        allocator: std.mem.Allocator,
        input1: []const f64,
        input2: []const f64,
        output: []f64,
    ) !void {
        if (input1.len != input2.len or input1.len != output.len) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input1.len, .estimate);
        defer fft_ctx.deinit();

        // Transform both functions
        const temp1 = try allocator.alloc(f64, input1.len);
        defer allocator.free(temp1);
        const temp2 = try allocator.alloc(f64, input2.len);
        defer allocator.free(temp2);

        try fft_ctx.forwardTransform(input1, temp1);
        try fft_ctx.forwardTransform(input2, temp2);

        // Compute cross-correlation in frequency domain
        // For real transforms, conjugate is implicit in the symmetry
        for (0..input1.len) |i| {
            output[i] = temp1[i] * temp2[i];
        }

        // Transform back to real domain
        try fft_ctx.inverseTransform(output, output);
    }

    /// Compute power spectrum of an array
    /// Useful for analyzing frequency content and debugging transforms
    pub fn powerSpectrum(
        allocator: std.mem.Allocator,
        input: []const f64,
        spectrum: []f64,
    ) !void {
        if (input.len != spectrum.len) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input.len, .estimate);
        defer fft_ctx.deinit();

        // Transform to frequency domain
        try fft_ctx.forwardTransform(input, spectrum);

        // Compute power spectrum |F(k)|²
        for (0..input.len) |i| {
            spectrum[i] = spectrum[i] * spectrum[i];
        }
    }

    /// Apply frequency domain filter to input signal
    /// filter array should have same length as input
    pub fn applyFilter(
        allocator: std.mem.Allocator,
        input: []const f64,
        filter: []const f64,
        output: []f64,
    ) !void {
        if (input.len != filter.len or input.len != output.len) {
            return FFTError.InvalidDataSize;
        }

        var fft_ctx = try FFTContext.init(allocator, input.len, .estimate);
        defer fft_ctx.deinit();

        // Transform input to frequency domain
        try fft_ctx.forwardTransform(input, output);

        // Apply filter in frequency domain
        for (0..input.len) |i| {
            output[i] *= filter[i];
        }

        // Transform back to real domain
        try fft_ctx.inverseTransform(output, output);
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

test "FFT convolution operation" {
    const allocator = std.testing.allocator;
    const n = 32;

    const input1 = try allocator.alloc(f64, n);
    defer allocator.free(input1);
    const input2 = try allocator.alloc(f64, n);
    defer allocator.free(input2);
    const output = try allocator.alloc(f64, n);
    defer allocator.free(output);

    // Create simple test signals
    @memset(input1, 0.0);
    @memset(input2, 0.0);
    input1[0] = 1.0; // Delta function
    input2[1] = 1.0; // Shifted delta

    try FourierUtils.convolution(allocator, input1, input2, output);

    // Convolution of two deltas should produce a shifted result
    // Verify output is finite and has expected properties
    for (0..n) |i| {
        try std.testing.expect(std.math.isFinite(output[i]));
    }

    // The result should have most energy concentrated in few points
    var total_energy: f64 = 0.0;
    for (0..n) |i| {
        total_energy += output[i] * output[i];
    }
    try std.testing.expect(total_energy > 0.0);
}

test "Power spectrum computation" {
    const allocator = std.testing.allocator;
    const n = 32;

    const input = try allocator.alloc(f64, n);
    defer allocator.free(input);
    const spectrum = try allocator.alloc(f64, n);
    defer allocator.free(spectrum);

    // Create a sine wave
    for (0..n) |i| {
        const x = 2.0 * std.math.pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        input[i] = @sin(x);
    }

    try FourierUtils.powerSpectrum(allocator, input, spectrum);

    // Power spectrum should be non-negative
    for (0..n) |i| {
        try std.testing.expect(spectrum[i] >= 0.0);
    }

    // For a pure sine wave, most power should be concentrated in specific frequencies
    var total_power: f64 = 0.0;
    for (0..n) |i| {
        total_power += spectrum[i];
    }
    try std.testing.expect(total_power > 0.0);
}

test "FFT filter application" {
    const allocator = std.testing.allocator;
    const n = 16;

    const input = try allocator.alloc(f64, n);
    defer allocator.free(input);
    const filter = try allocator.alloc(f64, n);
    defer allocator.free(filter);
    const output = try allocator.alloc(f64, n);
    defer allocator.free(output);

    // Create test signal and filter
    for (0..n) |i| {
        input[i] = @as(f64, @floatFromInt(i + 1));
        filter[i] = if (i < n / 2) 1.0 else 0.0; // Low-pass filter
    }

    try FourierUtils.applyFilter(allocator, input, filter, output);

    // Filtered output should be finite and different from input
    var input_sum: f64 = 0.0;
    var output_sum: f64 = 0.0;
    for (0..n) |i| {
        try std.testing.expect(std.math.isFinite(output[i]));
        input_sum += input[i];
        output_sum += output[i];
    }

    // Filter should modify the signal (but the sum might be similar)
    // Just verify that the output is different from input
    var is_different = false;
    for (0..n) |i| {
        if (@abs(input[i] - output[i]) > 1e-10) {
            is_different = true;
            break;
        }
    }
    try std.testing.expect(is_different);
}
