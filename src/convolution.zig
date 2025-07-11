const std = @import("std");
const math = std.math;
const testing = std.testing;
const dsbmv = @import("linear/dsbmv.zig");
const sb = @import("linear/sb.zig");

pub const WeightRule = enum {
    rectangular, // Uniform weights (all 1.0)
    simpson, // Simpson's rule weights
};

/// Integration weights vector for radially symmetric convolution kernels.
/// This should be used with a kernel vector whose first element is the center of the kernel,
/// and the rest are positive offsets from the center.
pub const RadialWeights = struct {
    data: []f64,

    const Self = @This();

    pub fn init(
        comptime rule: WeightRule,
        allocator: std.mem.Allocator,
        size: usize,
        spacing: f64,
    ) !RadialWeights {
        const data = try allocator.alloc(f64, size);

        if (rule == .simpson) {
            return try simpson(data, size, spacing);
        } else if (rule == .rectangular) {
            return try rectangular(data);
        } else {
            unreachable;
        }
    }

    fn simpson(
        data: []f64,
        n: usize,
        h: f64,
    ) !RadialWeights {
        const h_over_3: f64 = h / 3.0;

        // The first point is the center of the integration interval
        var coefficient: f64 = if (n % 2 == 0) 4.0 else 2.0;

        for (0..n) |i| {
            if (i == 0) {
                coefficient = if (n % 2 == 0) 4.0 else 2.0;
            } else if (i == n - 1) {
                coefficient = 1.0; // Last point always has coefficient 1
            } else {
                coefficient = if (coefficient == 4) 2.0 else 4.0;
            }
            data[i] = coefficient * h_over_3;
        }
        return RadialWeights{ .data = data };
    }

    fn rectangular(data: []f64) !RadialWeights {
        // Uniform weights (all 1.0) for rectangular rule
        @memset(data, 1.0);
        return RadialWeights{ .data = data };
    }

    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }
};

test "simpson_weights_2" {
    var weights = try RadialWeights.init(.simpson, testing.allocator, 2, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_3" {
    var weights = try RadialWeights.init(.simpson, testing.allocator, 3, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_4" {
    var weights = try RadialWeights.init(.simpson, testing.allocator, 4, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 4.0, 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_5" {
    var weights = try RadialWeights.init(.simpson, testing.allocator, 5, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 2.0, 4.0, 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_6" {
    var weights = try RadialWeights.init(.simpson, testing.allocator, 6, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 4.0, 2.0, 4.0, 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

/// Convolution kernel for symmetric band matrix operations
pub const Kernel = struct {
    matrix: sb.SymmetricBandMatrix,
    radius: usize,

    const Self = @This();

    /// Initialize a convolution kernel with integration weights applied
    /// kernel_values contains only half spectrum (center + positive offsets)
    /// weights are applied element-wise to the kernel values
    pub fn init(
        allocator: std.mem.Allocator,
        values: []const f64,
        signal_size: usize,
        weights: RadialWeights,
    ) !Self {
        const radius = values.len - 1;
        var matrix = try sb.SymmetricBandMatrix.init(allocator, signal_size, radius);
        matrix.clear();

        // Apply weights to kernel values and construct band matrix
        for (0..signal_size) |i| {
            for (0..@min(radius + 1, signal_size - i)) |offset| {
                const j = i + offset;
                const weighted_value = values[offset] * weights.data[offset];
                matrix.set(i, j, weighted_value);
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
        signal: []const f64,
        result: []f64,
    ) void {
        dsbmv.dsbmv(.U, self.matrix.n, self.matrix.k, 1.0, self.matrix.data, self.matrix.lda, signal, 1, 0.0, result, 1);
    }

    /// Deinitialize and free the kernel matrix
    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        self.matrix.deinit(allocator);
    }
};

test "convolve_triangular_kernel" {
    const allocator = testing.allocator;
    const n = 17;
    const h = 1.0 / @as(f64, @floatFromInt(n - 1));

    const weights = try RadialWeights.init(.simpson, allocator, n, h);
    defer weights.deinit(allocator);

    var kernel_values = try allocator.alloc(f64, n);
    defer allocator.free(kernel_values);
    for (0..n) |i| {
        kernel_values[i] = 1.0 - h * @as(f64, @floatFromInt(i));
    }
    try testing.expectApproxEqAbs(0.0, kernel_values[n - 1], 1e-6);

    const kernel = try Kernel.init(allocator, kernel_values, 32, weights);
    defer kernel.deinit(allocator);

    const signal = try allocator.alloc(f64, 32);
    defer allocator.free(signal);
    @memset(signal, 1.0);

    const result = try allocator.alloc(f64, 32);
    defer allocator.free(result);

    kernel.convolve(signal, result);
    try testing.expectApproxEqAbs(1.0, result[n - 2], 1e-6);
}
