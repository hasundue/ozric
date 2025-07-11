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

pub const Weights = struct {
    data: []f64,

    const Self = @This();

    pub fn init(
        comptime rule: WeightRule,
        allocator: Allocator,
        nodes: []const usize,
        h: f64,
    ) !Weights {
        const max_node = nodes[nodes.len - 1];
        const n = max_node + 1;
        const data = try allocator.alloc(f64, n);

        if (rule == .simpson) {
            return try initSimpsonWeights(data, h, nodes);
        } else if (rule == .rectangular) {
            return try initRectangularWeights(data);
        } else {
            unreachable;
        }
    }

    fn initSimpsonWeights(
        data: []f64,
        h: f64,
        nodes: []const usize,
    ) !Weights {
        const n = data.len;
        if (n < 3) return error.InsufficientPoints;
        if ((n - 1) % 2 != 0) return error.OddIntervals; // n-1 intervals must be even

        // Initialize all weights to zero
        @memset(data, 0.0);

        // Process each chunk directly without allocation, always starting from 0
        var prev_node: usize = 0;
        for (nodes) |node| {
            const start = prev_node;
            const end = node;
            const chunk_len = end - start + 1;

            if (chunk_len < 3) {
                prev_node = node;
                continue; // Skip chunks too small for Simpson's rule
            }

            // Check if chunk has even number of intervals (odd number of points)
            const intervals = chunk_len - 1;
            if (intervals % 2 != 0) {
                prev_node = node;
                continue; // Skip chunks with odd intervals
            }

            // Apply Simpson's rule: [1, 4, 2, 4, 2, ..., 4, 1] * h/3
            data[start] += h / 3.0;
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

    fn initRectangularWeights(data: []f64) !Weights {
        // Uniform weights (all 1.0) for rectangular rule
        @memset(data, 1.0);
        return Weights{ .data = data };
    }

    pub fn deinit(self: Self, allocator: Allocator) void {
        allocator.free(self.data);
    }
};

pub fn join(
    allocator: Allocator,
    a: Weights,
    b: Weights,
) !Weights {
    // concat two Weights at a shared discontinuity point
    // The last point of self and first point of other are at the same location (discontinuity)
    // So we add them together and create a new combined weight array
    const combined_len = a.data.len + b.data.len - 1;
    var combined_data = try allocator.alloc(f64, combined_len);

    // Copy first weight array
    @memcpy(combined_data[0..a.data.len], a.data);

    // Add the overlapping point (discontinuity gets double weight)
    combined_data[a.data.len - 1] += b.data[0];

    // Copy the rest of the second weight array
    @memcpy(combined_data[a.data.len..], b.data[1..]);

    return Weights{ .data = combined_data };
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
        const nodes = [_]usize{kernel.len - 1};
        const rect_weights = try Weights.init(.rectangular, allocator, &nodes, 0.1);
        defer rect_weights.deinit(allocator);
        var kernel_obj = try Kernel.init(allocator, kernel, size, rect_weights);
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
    const nodes = [_]usize{kernel.len - 1};
    const rect_weights = try Weights.init(.rectangular, allocator, &nodes, 0.1);
    defer rect_weights.deinit(allocator);
    var kernel_obj = try Kernel.init(allocator, &kernel, signal.len, rect_weights);
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
        const nodes1 = [_]usize{kernel1.len - 1};
        const rect_weights1 = try Weights.init(.rectangular, allocator, &nodes1, 0.1);
        defer rect_weights1.deinit(allocator);
        var kernel_obj = try Kernel.init(allocator, &kernel1, signal.len, rect_weights1);
        defer kernel_obj.deinit(allocator);
        kernel_obj.convolve(&signal, result);
        const expected_smooth = [_]f64{ 1.0, 2.0, 3.0, 4.0, 3.5 };
        try testing.expectEqualSlices(f64, expected_smooth[0..2], result[0..2]);
    }

    // Second convolution with identity kernel (should return original signal)
    {
        const nodes2 = [_]usize{kernel2.len - 1};
        const rect_weights2 = try Weights.init(.rectangular, allocator, &nodes2, 0.1);
        defer rect_weights2.deinit(allocator);
        var kernel_obj = try Kernel.init(allocator, &kernel2, signal.len, rect_weights2);
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

    const nodes = [_]usize{kernel.len - 1};
    const rect_weights = try Weights.init(.rectangular, allocator, &nodes, 0.1);
    defer rect_weights.deinit(allocator);
    var kernel_obj = try Kernel.init(allocator, &kernel, signal_size, rect_weights);
    defer kernel_obj.deinit(allocator);

    // Check that the kernel values are properly stored
    const test_positions = [_][2]usize{ .{ 0, 0 }, .{ 0, 1 }, .{ 1, 0 }, .{ 1, 1 }, .{ 1, 2 } };
    const expected_values = [_]f64{ 0.5, 0.25, 0.25, 0.5, 0.25 };
    for (test_positions, expected_values) |pos, expected| {
        try eqa(expected, kernel_obj.matrix.get(pos[0], pos[1]), 1e-10);
    }
}

test "Simpson weights" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Test Simpson weights for 5 points (4 intervals)
    const nodes = [_]usize{4};
    var weights = try Weights.init(.simpson, allocator, &nodes, 0.1);
    defer weights.deinit(allocator);

    // Expected pattern: [1, 4, 2, 4, 1] * (0.1/3)
    const expected = [_]f64{ 0.1 / 3.0, 4 * 0.1 / 3.0, 2 * 0.1 / 3.0, 4 * 0.1 / 3.0, 0.1 / 3.0 };
    try testing.expectEqual(5, weights.data.len);
    for (expected, weights.data) |exp, actual| {
        try eqa(exp, actual, 1e-10);
    }
}

test "Weights join" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Create two weight arrays that should be concatenated at a discontinuity
    const nodes1 = [_]usize{2};
    const nodes2 = [_]usize{2};
    var weights1 = try Weights.init(.simpson, allocator, &nodes1, 0.1);
    var weights2 = try Weights.init(.simpson, allocator, &nodes2, 0.1);

    // Join them
    var combined = try join(allocator, weights1, weights2);
    defer combined.deinit(allocator);
    defer weights1.deinit(allocator);
    defer weights2.deinit(allocator);

    // Should have 5 elements (3 + 3 - 1)
    try testing.expectEqual(5, combined.data.len);

    // Middle element should be sum of overlapping weights
    const expected_middle = weights1.data[2] + weights2.data[0];
    try eqa(expected_middle, combined.data[2], 1e-10);
}

test "Kernel with weights" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Create kernel values and weights (Simpson needs at least 3 points)
    const kernel_values = [_]f64{ 1.0, 0.5, 0.25 };
    const nodes = [_]usize{2};
    var weights = try Weights.init(.simpson, allocator, &nodes, 0.1);
    defer weights.deinit(allocator);

    // Create kernel with weights applied
    var kernel = try Kernel.init(allocator, &kernel_values, 4, weights);
    defer kernel.deinit(allocator);

    // Check that weights were applied to kernel values
    // Simpson weights for 3 points: [1, 4, 1] * (0.1/3) = [0.1/3, 0.4/3, 0.1/3]
    const expected_center = kernel_values[0] * weights.data[0];
    const expected_offset1 = kernel_values[1] * weights.data[1];
    const expected_offset2 = kernel_values[2] * weights.data[2];

    try eqa(expected_center, kernel.matrix.get(0, 0), 1e-10);
    try eqa(expected_offset1, kernel.matrix.get(0, 1), 1e-10);
    try eqa(expected_offset2, kernel.matrix.get(0, 2), 1e-10);
}

test "Simpson weights with kinks" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Test Simpson weights with a kink at index 2 (separating into two chunks)
    // Array: [0, 1, 2, 3, 4] with nodes at [2, 4] (0 is implicit)
    // Chunk 1: [0, 2] -> Simpson: [1, 4, 1] * (0.1/3)
    // Chunk 2: [2, 4] -> Simpson: [1, 4, 1] * (0.1/3)
    // Result: [0.1/3, 4*0.1/3, 2*0.1/3, 4*0.1/3, 0.1/3]
    //         (index 2 gets weight from both chunks)
    const nodes = [_]usize{ 2, 4 };
    var weights = try Weights.init(.simpson, allocator, &nodes, 0.1);
    defer weights.deinit(allocator);

    const expected = [_]f64{
        0.1 / 3.0, // start of chunk 1
        4.0 * 0.1 / 3.0, // middle of chunk 1
        2.0 * 0.1 / 3.0, // end of chunk 1 + start of chunk 2
        4.0 * 0.1 / 3.0, // middle of chunk 2
        0.1 / 3.0, // end of chunk 2
    };
    try testing.expectEqual(5, weights.data.len);
    for (expected, weights.data) |exp, actual| {
        try eqa(exp, actual, 1e-10);
    }
}

test "Simpson weights with multiple kinks" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Test with 7 points and nodes at indices 2, 4, 6 (0 is implicit)
    // Array: [0, 1, 2, 3, 4, 5, 6] with nodes at [2, 4, 6] (0 is implicit)
    // Chunk 1: [0, 2] -> Simpson: [1, 4, 1] * (0.1/3)
    // Chunk 2: [2, 4] -> Simpson: [1, 4, 1] * (0.1/3)
    // Chunk 3: [4, 6] -> Simpson: [1, 4, 1] * (0.1/3)
    const nodes = [_]usize{ 2, 4, 6 };
    var weights = try Weights.init(.simpson, allocator, &nodes, 0.1);
    defer weights.deinit(allocator);

    const expected = [_]f64{
        0.1 / 3.0, // start of chunk 1
        4.0 * 0.1 / 3.0, // middle of chunk 1
        2.0 * 0.1 / 3.0, // end of chunk 1 + start of chunk 2
        4.0 * 0.1 / 3.0, // middle of chunk 2
        2.0 * 0.1 / 3.0, // end of chunk 2 + start of chunk 3
        4.0 * 0.1 / 3.0, // middle of chunk 3
        0.1 / 3.0, // end of chunk 3
    };
    try testing.expectEqual(7, weights.data.len);
    for (expected, weights.data) |exp, actual| {
        try eqa(exp, actual, 1e-10);
    }
}
