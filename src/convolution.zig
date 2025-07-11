const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;
const dsbmv = @import("linear/dsbmv.zig");
const sb = @import("linear/sb.zig");

pub const WeightRule = enum {
    rectangular, // Uniform weights (all 1.0)
    simpson, // Simpson's rule weights
};

/// Integration weights for radially symmetric convolution kernels.
/// This should be used with a kernel vector whose first element is the center of the kernel,
/// and the rest are positive offsets from the center.
pub const Weights = struct {
    data: []f64,

    const Self = @This();

    pub fn init(
        comptime rule: WeightRule,
        allocator: Allocator,
        /// Indices of kinks and the edge of the kernel
        nodes: []const usize,
        /// Spacing of the kernel vector
        h: f64,
    ) !Weights {
        const max_node = nodes[nodes.len - 1];
        const n = max_node + 1;
        const data = try allocator.alloc(f64, n);

        if (rule == .simpson) {
            return try initSimpson(data, h, nodes);
        } else if (rule == .rectangular) {
            return try initRectangular(data);
        } else {
            unreachable;
        }
    }

    fn initSimpson(
        data: []f64,
        h: f64,
        nodes: []const usize,
    ) !Weights {
        // Initialize all weights to zero
        @memset(data, 0.0);

        // Process each chunk directly without allocation, always starting from 0
        var prev_node: usize = 0;
        for (nodes) |node| {
            const start = prev_node;
            const end = node;
            const chunk_len = end - start + 1;

            if (start == 0) {
                // The first point is the center of the actual integration interval
                const coefficient: f64 = if (chunk_len % 2 == 0) 4.0 else 2.0;
                data[start] += coefficient * h / 3.0;
            } else {
                data[start] += h / 3.0;
            }
            data[end] += h / 3.0;

            for (start + 1..end) |i| {
                const offset_in_chunk = i - start;
                const coefficient: f64 = if (offset_in_chunk % 2 == 1) 4.0 else 2.0;
                data[i] += coefficient * h / 3.0;
            }

            prev_node = node;
        }

        return Weights{ .data = data };
    }

    fn initRectangular(data: []f64) !Weights {
        // Uniform weights (all 1.0) for rectangular rule
        @memset(data, 1.0);
        return Weights{ .data = data };
    }

    pub fn deinit(self: Self, allocator: Allocator) void {
        allocator.free(self.data);
    }
};

test "simpson_weights_single_node_1" {
    var weights = try Weights.init(.simpson, testing.allocator, &[_]usize{1}, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_single_node_2" {
    var weights = try Weights.init(.simpson, testing.allocator, &[_]usize{2}, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_multiple_nodes_1" {
    var weights = try Weights.init(.simpson, testing.allocator, &[_]usize{ 1, 3 }, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 4.0, 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_multiple_nodes_2" {
    var weights = try Weights.init(.simpson, testing.allocator, &[_]usize{ 2, 4 }, 3.0);
    defer weights.deinit(testing.allocator);
    const expected = [_]f64{ 2.0, 4.0, 2.0, 4.0, 1.0 };
    for (expected, weights.data) |exp, actual| try testing.expectApproxEqAbs(exp, actual, 1e-10);
}

test "simpson_weights_multiple_nodes_3" {
    var weights = try Weights.init(.simpson, testing.allocator, &[_]usize{ 1, 5 }, 3.0);
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
        allocator: Allocator,
        kernel_values: []const f64,
        signal_size: usize,
        weights: Weights,
    ) !Self {
        const radius = kernel_values.len - 1;
        var matrix = try sb.SymmetricBandMatrix.init(allocator, signal_size, radius);
        matrix.clear();

        // Apply weights to kernel values and construct band matrix
        for (0..signal_size) |i| {
            for (0..@min(radius + 1, signal_size - i)) |offset| {
                const j = i + offset;
                const weighted_value = kernel_values[offset] * weights.data[offset];
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
    pub fn deinit(self: Self, allocator: Allocator) void {
        self.matrix.deinit(allocator);
    }
};

test "convolve_triangular_kernel" {
    const allocator = testing.allocator;
    const n = 17;
    const h = 1.0 / @as(f64, @floatFromInt(n - 1));

    const weights = try Weights.init(.simpson, allocator, &.{n - 1}, h);
    defer weights.deinit(allocator);

    var kernel_values = try allocator.alloc(f64, n);
    defer allocator.free(kernel_values);
    for (0..n) |i| {
        kernel_values[i] = 1.0 - h * @as(f64, @floatFromInt(i));
    }

    const kernel = try Kernel.init(allocator, kernel_values, 2 * n, weights);
    defer kernel.deinit(allocator);

    const values = try allocator.alloc(f64, 2 * n);
    defer allocator.free(values);
    @memset(values, 1.0);

    const result = try allocator.alloc(f64, 2 * n);
    defer allocator.free(result);

    kernel.convolve(values, result);
    try testing.expectApproxEqAbs(1.0, result[n - 1], 1e-6);
}
