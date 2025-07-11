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

test "rectangular convolution simpson vs rectangular" {
    const allocator = testing.allocator;

    // Test convolution of two rectangular functions
    // Both signal and kernel are smooth, so Simpson's rule should be more accurate
    const dt = 0.1;
    const n = 65; // Increased grid resolution for better accuracy testing

    // Create rectangular signal: f(t) = 1 for |t| ≤ 2.0, 0 otherwise
    // Domain extends from -3.2 to 3.2 with dt=0.1 for n=65
    // Stretch signal to better utilize the grid
    var signal = try allocator.alloc(f64, n);
    defer allocator.free(signal);
    for (0..n) |i| {
        const t = -3.2 + @as(f64, @floatFromInt(i)) * dt;
        signal[i] = if (@abs(t) <= 2.0) 1.0 else 0.0;
    }

    // Create rectangular kernel: g(t) = 1 for |t| ≤ 1.6, 0 otherwise
    // Half-spectrum representation: [g(0), g(dt), g(2*dt), ...]
    // Stretch kernel to better utilize the grid resolution
    // Need odd number of points for Simpson's rule (even number of intervals)
    const kernel_half_width = 16; // 1.6/0.1 = 16 points, giving 17 total points
    var rect_kernel = try allocator.alloc(f64, kernel_half_width + 1);
    defer allocator.free(rect_kernel);
    for (0..rect_kernel.len) |i| {
        const t = @as(f64, @floatFromInt(i)) * dt;
        rect_kernel[i] = if (t <= 1.6) 1.0 else 0.0;
    }

    const result_rect = try allocator.alloc(f64, n);
    defer allocator.free(result_rect);
    const result_simpson = try allocator.alloc(f64, n);
    defer allocator.free(result_simpson);

    // Rectangular rule
    {
        const nodes = [_]usize{rect_kernel.len - 1};
        const weights = try Weights.init(.rectangular, allocator, &nodes, dt);
        defer weights.deinit(allocator);
        var kernel = try Kernel.init(allocator, rect_kernel, n, weights);
        defer kernel.deinit(allocator);
        kernel.convolve(signal, result_rect);
    }

    // Simpson's rule
    {
        const nodes = [_]usize{rect_kernel.len - 1};
        const weights = try Weights.init(.simpson, allocator, &nodes, dt);
        defer weights.deinit(allocator);
        var kernel = try Kernel.init(allocator, rect_kernel, n, weights);
        defer kernel.deinit(allocator);
        kernel.convolve(signal, result_simpson);
    }

    // Analytical solution: convolution of two rectangles is a triangular function
    // For rect(-2.0,2.0) * rect(-1.6,1.6), the result is triangular with max at t=0
    // and width 3.6 (sum of half-widths: 2.0 + 1.6 = 3.6)
    var analytical = try allocator.alloc(f64, n);
    defer allocator.free(analytical);
    for (0..n) |i| {
        const t = -3.2 + @as(f64, @floatFromInt(i)) * dt;
        if (@abs(t) <= 3.6) {
            analytical[i] = 3.6 - @abs(t); // triangular function
        } else {
            analytical[i] = 0.0;
        }
    }

    // Print comparison
    std.debug.print("\nRectangular convolution: rect(-2.0,2.0) * rect(-1.6,1.6)\n", .{});
    std.debug.print("t\tAnalytical\tRect\tSimpson\tRect_err\tSimp_err\n", .{});
    var max_rect_err: f64 = 0.0;
    var max_simp_err: f64 = 0.0;
    for (0..n) |i| {
        const t = -3.2 + @as(f64, @floatFromInt(i)) * dt;
        const expected = analytical[i];
        const rect_err = @abs(result_rect[i] - expected);
        const simp_err = @abs(result_simpson[i] - expected);
        max_rect_err = @max(max_rect_err, rect_err);
        max_simp_err = @max(max_simp_err, simp_err);
        std.debug.print("{d:.1}\t{d:.3}\t\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.3}\n", .{ t, expected, result_rect[i], result_simpson[i], rect_err, simp_err });
    }

    std.debug.print("Max errors: Rectangular={d:.3}, Simpson={d:.3}\n", .{ max_rect_err, max_simp_err });

    // Simpson's rule should be more accurate for smooth functions
    try testing.expect(max_simp_err < max_rect_err);

    // Verify the dramatic improvement in accuracy
    try testing.expect(max_rect_err > 3.0); // Rectangular rule has large errors
    try testing.expect(max_simp_err < 0.5); // Simpson's rule is much more accurate
}

test "simpson accuracy vs grid resolution" {
    const allocator = testing.allocator;

    // Test Simpson's rule accuracy for rect * rect convolution at different grid resolutions
    const resolutions = [_]usize{ 17, 33, 65, 129 }; // Powers of 2 + 1 for Simpson's rule

    std.debug.print("\nSimpson's rule accuracy vs grid resolution:\n", .{});
    std.debug.print("N\tdt\tMax_Error\tL2_Error\tConv_Rate\n", .{});

    var prev_max_error: f64 = 0.0;
    var prev_dt: f64 = 0.0;

    for (resolutions, 0..) |n, i| {
        // Domain: -3.2 to 3.2, adjust dt based on n
        const dt = 6.4 / @as(f64, @floatFromInt(n - 1));

        // Create rectangular signal: f(t) = 1 for |t| ≤ 2.0, 0 otherwise
        var signal = try allocator.alloc(f64, n);
        defer allocator.free(signal);
        for (0..n) |j| {
            const t = -3.2 + @as(f64, @floatFromInt(j)) * dt;
            signal[j] = if (@abs(t) <= 2.0) 1.0 else 0.0;
        }

        // Create rectangular kernel: g(t) = 1 for |t| ≤ 1.6, 0 otherwise
        // Scale kernel size based on dt to maintain same physical size
        const kernel_half_width = @as(usize, @intFromFloat(1.6 / dt));
        // Ensure odd number of points for Simpson's rule
        const actual_kernel_size = if (kernel_half_width % 2 == 0) kernel_half_width + 1 else kernel_half_width;

        var rect_kernel = try allocator.alloc(f64, actual_kernel_size);
        defer allocator.free(rect_kernel);
        for (0..rect_kernel.len) |j| {
            const t = @as(f64, @floatFromInt(j)) * dt;
            rect_kernel[j] = if (t <= 1.6) 1.0 else 0.0;
        }

        const result_simpson = try allocator.alloc(f64, n);
        defer allocator.free(result_simpson);

        // Simpson's rule convolution
        {
            const nodes = [_]usize{rect_kernel.len - 1};
            const weights = try Weights.init(.simpson, allocator, &nodes, dt);
            defer weights.deinit(allocator);
            var kernel = try Kernel.init(allocator, rect_kernel, n, weights);
            defer kernel.deinit(allocator);
            kernel.convolve(signal, result_simpson);
        }

        // Analytical solution: triangular function
        var analytical = try allocator.alloc(f64, n);
        defer allocator.free(analytical);
        for (0..n) |j| {
            const t = -3.2 + @as(f64, @floatFromInt(j)) * dt;
            if (@abs(t) <= 3.6) {
                analytical[j] = 3.6 - @abs(t);
            } else {
                analytical[j] = 0.0;
            }
        }

        // Compute errors
        var max_error: f64 = 0.0;
        var l2_error: f64 = 0.0;
        for (0..n) |j| {
            const err = @abs(result_simpson[j] - analytical[j]);
            max_error = @max(max_error, err);
            l2_error += err * err;
        }
        l2_error = @sqrt(l2_error / @as(f64, @floatFromInt(n)));

        // Compute convergence rate (should be ~4th order for Simpson's rule)
        var convergence_rate: f64 = 0.0;
        if (i > 0 and prev_max_error > 0.0) {
            const dt_ratio = prev_dt / dt;
            convergence_rate = @log(prev_max_error / max_error) / @log(dt_ratio);
        }

        std.debug.print("{d}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.2}\n", .{ n, dt, max_error, l2_error, convergence_rate });

        prev_max_error = max_error;
        prev_dt = dt;
    }

    std.debug.print("Note: Convergence rate should approach 4.0 for Simpson's rule\n", .{});
}

test "simpson gaussian convolution convergence" {
    const allocator = testing.allocator;

    // Test Simpson's rule accuracy for Gaussian * Gaussian convolution at different grid resolutions
    // Gaussian convolution is perfectly smooth and should show 4th-order convergence
    const resolutions = [_]usize{ 17, 33, 65, 129 }; // Powers of 2 + 1 for Simpson's rule

    std.debug.print("\nSimpson's rule for Gaussian * Gaussian convolution:\n", .{});
    std.debug.print("N\tdt\tMax_Error\tL2_Error\tConv_Rate\n", .{});

    var prev_max_error: f64 = 0.0;
    var prev_dt: f64 = 0.0;

    for (resolutions, 0..) |n, i| {
        // Domain: -4.0 to 4.0, adjust dt based on n
        const dt = 8.0 / @as(f64, @floatFromInt(n - 1));

        // Create Gaussian signal: f(t) = exp(-t²/2σ²), σ = 1.0
        var signal = try allocator.alloc(f64, n);
        defer allocator.free(signal);
        for (0..n) |j| {
            const t = -4.0 + @as(f64, @floatFromInt(j)) * dt;
            signal[j] = @exp(-0.5 * t * t); // σ = 1.0
        }

        // Create Gaussian kernel: g(t) = exp(-t²/2σ²), σ = 0.8
        // Scale kernel size to cover ~3σ = 2.4
        const kernel_half_width = @as(usize, @intFromFloat(2.4 / dt));
        // Ensure odd number of points for Simpson's rule
        const actual_kernel_size = if (kernel_half_width % 2 == 0) kernel_half_width + 1 else kernel_half_width;

        var gauss_kernel = try allocator.alloc(f64, actual_kernel_size);
        defer allocator.free(gauss_kernel);
        for (0..gauss_kernel.len) |j| {
            const t = @as(f64, @floatFromInt(j)) * dt;
            gauss_kernel[j] = @exp(-0.5 * t * t / (0.8 * 0.8)); // σ = 0.8
        }

        const result_simpson = try allocator.alloc(f64, n);
        defer allocator.free(result_simpson);

        // Simpson's rule convolution
        {
            const nodes = [_]usize{gauss_kernel.len - 1};
            const weights = try Weights.init(.simpson, allocator, &nodes, dt);
            defer weights.deinit(allocator);
            var kernel = try Kernel.init(allocator, gauss_kernel, n, weights);
            defer kernel.deinit(allocator);
            kernel.convolve(signal, result_simpson);
        }

        // Analytical solution: Gaussian * Gaussian = Gaussian with σ_result = √(σ₁² + σ₂²)
        // For σ₁ = 1.0, σ₂ = 0.8: σ_result = √(1.0 + 0.64) = √1.64 ≈ 1.281
        const sigma_result = @sqrt(1.0 + 0.8 * 0.8);
        const amplitude = @sqrt(2.0 * math.pi) * 0.8; // Normalization factor

        var analytical = try allocator.alloc(f64, n);
        defer allocator.free(analytical);
        for (0..n) |j| {
            const t = -4.0 + @as(f64, @floatFromInt(j)) * dt;
            analytical[j] = amplitude * @exp(-0.5 * t * t / (sigma_result * sigma_result));
        }

        // Compute errors
        var max_error: f64 = 0.0;
        var l2_error: f64 = 0.0;
        for (0..n) |j| {
            const err = @abs(result_simpson[j] - analytical[j]);
            max_error = @max(max_error, err);
            l2_error += err * err;
        }
        l2_error = @sqrt(l2_error / @as(f64, @floatFromInt(n)));

        // Compute convergence rate (should approach 4.0 for smooth functions)
        var convergence_rate: f64 = 0.0;
        if (i > 0 and prev_max_error > 0.0) {
            const dt_ratio = prev_dt / dt;
            convergence_rate = @log(prev_max_error / max_error) / @log(dt_ratio);
        }

        std.debug.print("{d}\t{d:.3}\t{d:.3}\t{d:.3}\t{d:.2}\n", .{ n, dt, max_error, l2_error, convergence_rate });

        prev_max_error = max_error;
        prev_dt = dt;
    }

    std.debug.print("Note: Convergence rate should approach 4.0 for smooth Gaussian functions\n", .{});
}
