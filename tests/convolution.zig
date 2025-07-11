const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;

const root = @import("convolution");
const convolution = root.convolution;
const Weights = convolution.Weights;
const Kernel = convolution.Kernel;

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

test "delta convolution simpson vs rectangular" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Test polynomial: f(t) = t² + 2t + 1 = (t+1)²
    // Domain: t ∈ [-1, 1], discretized with dt = 0.2
    const dt = 0.2;
    const n = 11; // -1.0 to 1.0 with step 0.2

    // Create polynomial signal f(t) = (t+1)²
    var signal = try allocator.alloc(f64, n);
    defer allocator.free(signal);
    for (0..n) |i| {
        const t = -1.0 + @as(f64, @floatFromInt(i)) * dt;
        signal[i] = (t + 1.0) * (t + 1.0); // (t+1)²
    }

    // Create discrete delta function at center
    const delta_kernel = [_]f64{ 1.0, 0.0, 0.0, 0.0, 0.0 }; // Delta at t=0

    const result_rect = try allocator.alloc(f64, n);
    defer allocator.free(result_rect);
    const result_simpson = try allocator.alloc(f64, n);
    defer allocator.free(result_simpson);

    // Rectangular rule
    {
        const nodes = [_]usize{delta_kernel.len - 1};
        const weights = try Weights.init(.rectangular, allocator, &nodes, dt);
        defer weights.deinit(allocator);
        var kernel = try Kernel.init(allocator, &delta_kernel, n, weights);
        defer kernel.deinit(allocator);
        kernel.convolve(signal, result_rect);
    }

    // Simpson's rule
    {
        const nodes = [_]usize{delta_kernel.len - 1};
        const weights = try Weights.init(.simpson, allocator, &nodes, dt);
        defer weights.deinit(allocator);
        var kernel = try Kernel.init(allocator, &delta_kernel, n, weights);
        defer kernel.deinit(allocator);
        kernel.convolve(signal, result_simpson);
    }

    // Print comparison
    std.debug.print("\nDelta convolution: f(t) = (t+1)², (f * δ)(t) = f(t)\n", .{});
    std.debug.print("t\tf(t)\tRect\tSimpson\tRect_err\tSimp_err\n", .{});
    for (0..n) |i| {
        const t = -1.0 + @as(f64, @floatFromInt(i)) * dt;
        const expected = signal[i];
        const rect_err = @abs(result_rect[i] - expected);
        const simp_err = @abs(result_simpson[i] - expected);
        std.debug.print("{d:.1}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\n", .{ t, expected, result_rect[i], result_simpson[i], rect_err, simp_err });
    }

    // Verify rectangular rule gives exact result for delta convolution
    for (signal, result_rect) |expected, actual| {
        try eqa(expected, actual, 1e-10);
    }
}
