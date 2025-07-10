const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const dsbmv = @import("linear/dsbmv.zig");
const sb = @import("linear/sb.zig");

pub const WeightRule = enum {
    // .Trapezoidal, // Trapezoidal rule weights (not implemented)
    simpson, // Simpson's rule weights
};

pub const Weights = struct {
    data: []f64,

    const Self = @This();

    pub fn init(
        comptime rule: WeightRule,
        allocator: Allocator,
        n: usize,
        h: f64,
    ) !Weights {
        var data = try allocator.alloc(f64, n);

        if (rule == .simpson) {
            return try initSimpsonWeights(allocator, &data, h);
        } else {
            unreachable;
        }
    }

    fn initSimpsonWeights(
        allocator: Allocator,
        data: *[]f64,
        h: f64,
    ) !Weights {}

    pub fn concat(
        self: Self,
        other: Weights,
    ) !Weights {
        // concat two Weights at a shared discontinuity point
    }
};

/// Convolution kernel for symmetric band matrix operations
pub const Kernel = struct {
    matrix: sb.SymmetricBandMatrix,
    radius: usize,

    const Self = @This();

    /// Initialize a convolution kernel and create its symmetric band matrix
    /// kernel_values contains only half spectrum (center + positive offsets)
    pub fn init(
        allocator: Allocator,
        kernel_values: []const f64,
        signal_size: usize,
    ) !Self {
        const radius = kernel_values.len - 1;
        var matrix = try sb.SymmetricBandMatrix.init(allocator, signal_size, radius);
        matrix.clear();

        // Construct band matrix from symmetric kernel
        for (0..signal_size) |i| {
            for (0..@min(radius + 1, signal_size - i)) |offset| {
                const j = i + offset;
                if (j < signal_size and offset < kernel_values.len) {
                    matrix.set(i, j, kernel_values[offset]);
                }
            }
        }

        return Self{
            .matrix = matrix,
            .radius = radius,
        };
    }

    /// Perform convolution with the kernel
    pub fn convolve(
        self: *const Self,
        weights: Weights,
        signal: []const f64,
        result: []f64,
    ) void {
        // TODO: do result = self x weights here
        dsbmv.dsbmv(.U, self.matrix.n, self.matrix.k, 1.0, self.matrix.data, self.matrix.lda, signal, 1, 0.0, result, 1);
    }

    /// Deinitialize and free the kernel matrix
    pub fn deinit(self: *Self, allocator: Allocator) void {
        self.matrix.deinit(allocator);
    }
};

/// Benchmarking function for convolution using dsbmv
pub fn benchmark_convolution(allocator: Allocator) !void {
    const sizes = [_]usize{ 100, 500, 1000, 5000 };
    const kernel_radius = 10;

    for (sizes) |size| {
        var signal = try allocator.alloc(f64, size);
        defer allocator.free(signal);

        const result = try allocator.alloc(f64, size);
        defer allocator.free(result);

        var kernel = try allocator.alloc(f64, kernel_radius + 1);
        defer allocator.free(kernel);

        // Generation of Gaussian kernel (half spectrum)
        for (0..kernel_radius + 1) |i| {
            const x = @as(f64, @floatFromInt(i));
            kernel[i] = math.exp(-0.5 * x * x);
        }

        // Generation of test data
        for (0..size) |i| {
            signal[i] = @sin(@as(f64, @floatFromInt(i)) * 0.1);
        }

        // Create kernel and execute benchmark
        var kernel_obj = try Kernel.init(allocator, kernel, size);
        defer kernel_obj.deinit(allocator);

        const start_time = std.time.nanoTimestamp();

        kernel_obj.convolve(signal, result);

        const end_time = std.time.nanoTimestamp();
        const elapsed_ns = @as(f64, @floatFromInt(end_time - start_time));
        const elapsed_ms = elapsed_ns / 1_000_000.0;

        std.debug.print("Size: {}, Time: {d:.2} ms\n", .{ size, elapsed_ms });
    }
}

test "convolution test" {
    const allocator = testing.allocator;

    // Simple convolution test with known expected results
    const kernel = [_]f64{ 0.5, 0.25 }; // smoothing kernel (half spectrum: center + positive offset)
    const signal = [_]f64{ 1, 2, 3, 4, 5 };
    const result = try allocator.alloc(f64, signal.len);
    defer allocator.free(result);

    // Create kernel
    var kernel_obj = try Kernel.init(allocator, &kernel, signal.len);
    defer kernel_obj.deinit(allocator);

    kernel_obj.convolve(&signal, result);

    // Expected results for smoothing kernel [0.5, 0.25] (symmetric: [0.25, 0.5, 0.25])
    // Manual calculation considering boundary effects with symmetric band matrix:
    // The kernel [0.5, 0.25] represents center value 0.5 and offset value 0.25
    // This creates a symmetric kernel where A[i,i] = 0.5 and A[i,i+1] = A[i+1,i] = 0.25

    // Verify key results based on actual symmetric band matrix convolution behavior
    const expected = [_]f64{ 1.0, 2.0, 3.0, 4.0, 3.5 };
    try testing.expectEqualSlices(f64, &expected, result);

    // Verify that the result is reasonable (all values positive and bounded)
    for (result) |val| {
        try testing.expect(val >= 0.0 and val <= 10.0);
    }
}

test "kernel comparison" {
    const allocator = testing.allocator;

    // Test different kernels
    const kernel1 = [_]f64{ 0.5, 0.25 }; // smoothing kernel (half spectrum)
    const kernel2 = [_]f64{ 1.0, 0.0 }; // identity kernel (half spectrum)
    const signal = [_]f64{ 1, 2, 3, 4, 5 };
    const result = try allocator.alloc(f64, signal.len);
    defer allocator.free(result);

    // First convolution with smoothing kernel
    {
        var kernel_obj = try Kernel.init(allocator, &kernel1, signal.len);
        defer kernel_obj.deinit(allocator);
        kernel_obj.convolve(&signal, result);
        const expected_smooth = [_]f64{ 1.0, 2.0, 3.0, 4.0, 3.5 };
        try testing.expectEqualSlices(f64, expected_smooth[0..2], result[0..2]);
    }

    // Second convolution with identity kernel (should return original signal)
    {
        var kernel_obj = try Kernel.init(allocator, &kernel2, signal.len);
        defer kernel_obj.deinit(allocator);
        kernel_obj.convolve(&signal, result);
        try testing.expectEqualSlices(f64, &signal, result);
    }
}

test "convolution matrix creation" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    const kernel = [_]f64{ 0.5, 0.25 };
    const signal_size = 5;

    var kernel_obj = try Kernel.init(allocator, &kernel, signal_size);
    defer kernel_obj.deinit(allocator);

    // Check that the kernel values are properly stored
    const test_positions = [_][2]usize{ .{ 0, 0 }, .{ 0, 1 }, .{ 1, 0 }, .{ 1, 1 }, .{ 1, 2 } };
    const expected_values = [_]f64{ 0.5, 0.25, 0.25, 0.5, 0.25 };
    for (test_positions, expected_values) |pos, expected| {
        try eqa(expected, kernel_obj.matrix.get(pos[0], pos[1]), 1e-10);
    }
}
