const std = @import("std");
const math = std.math;
const t = std.testing;
const conv = @import("convolution.zig");

const Grid = @import("grid.zig").Grid;

/// Number of grid points per diameter for optimal discretization
const GRID_PPD: usize = 16;

/// Weight function with its grid node indices
pub const WeightFunction = struct {
    ptr: *const fn (HardSphereDFT, f64) f64,

    /// Grid node distances relative to the hard sphere diameter
    /// Array ends at the cutoff radius, with kinks as intermediate points
    nodes: []const usize,
};

pub const HardSphereDFT = struct {
    /// Hard-sphere diameter
    diameter: f64,

    /// Pre-calculated optimal grid spacing that aligns diameter with a grid point
    resolution: f64,

    const Self = @This();

    pub fn init(diameter: f64) Self {
        const resolution = diameter / @as(f64, @floatFromInt(GRID_PPD));
        return Self{
            .diameter = diameter,
            .resolution = resolution,
        };
    }

    /// Get weight functions with precalculated grid node indices
    pub fn getWeightFns(self: Self) [3]WeightFunction {
        _ = self; // Weight functions are normalized, diameter doesn't affect node indices
        return [_]WeightFunction{
            .{ .ptr = weightFn0, .nodes = &[_]usize{1} }, // 0 to 1σ (0 is implicit)
            .{ .ptr = weightFn1, .nodes = &[_]usize{ 1, 2 } }, // 0 to 2σ with kink at 1σ (0 is implicit)
            .{ .ptr = weightFn2, .nodes = &[_]usize{1} }, // 0 to 1σ (0 is implicit)
        };
    }

    pub fn weightFn0(self: Self, r: f64) f64 {
        if (r > self.diameter) return 0.0;

        const sigma3 = math.pow(f64, self.diameter, 3);
        const pi = math.pi;

        return 3 / (4 * pi * sigma3);
    }

    pub fn weightFn1(hs: HardSphereDFT, r: f64) f64 {
        if (r < hs.diameter) {
            return hs.weightFn1core(r);
        } else {
            return hs.weightFn1tail(r);
        }
    }

    // The first order weight function inside the hard sphere
    pub fn weightFn1core(self: Self, r: f64) f64 {
        if (r > self.diameter) return 0.0;

        const x = r / self.diameter;
        const x2 = math.pow(f64, x, 2);

        const a_0 = 0.90724;
        const a_1 = -1.23717;
        const a_2 = 0.21616;

        return a_0 + a_1 * x + a_2 * x2;
    }

    /// The second order weight function outside the hard sphere
    pub fn weightFn1tail(self: Self, r: f64) f64 {
        if (r < self.diameter) return 0.0;

        const x = r / self.diameter;
        const x2 = math.pow(f64, x, 2);
        const x3 = math.pow(f64, x, 3);

        const c = -0.10244;
        const b_0 = 35.134;
        const alpha = 4.934;
        const b_1 = -98.684;
        const beta_1 = 3.5621;
        const b_2 = 92.693;
        const beta_2 = 12.0;
        const b_3 = -29.257;

        const first = c * @exp(-beta_1 * (x - 1.0)) * @sin(alpha * (x - 1.0));
        const second = @exp(-beta_2 * (x - 1.0)) * (b_0 + b_1 * x + b_2 * x2 + b_3 * x3);

        return first + second;
    }

    pub fn weightFn2(self: Self, r: f64) f64 {
        if (r >= self.diameter) return 0.0;

        const x = r / self.diameter;
        const x2 = x * x;
        const pi = math.pi;

        return (15.0 / 8.0 / pi / self.diameter) * (1 - 3 * x + 3 * x2);
    }
};

test "HardSphereDFT init" {
    const hs = HardSphereDFT.init(1.0);
    try t.expectEqual(@as(f64, 1.0), hs.diameter);
    try t.expectEqual(@as(f64, 1.0 / 16.0), hs.resolution);
}

test "HardSphereDFT optimal spacing calculation" {
    const hs = HardSphereDFT.init(1.0);

    // With diameter=1.0 and GRID_PPD=16, resolution should be 1.0/16
    try t.expectEqual(@as(f64, 1.0 / 16.0), hs.resolution);

    // Verify that the hard sphere diameter falls exactly on a grid point
    // At resolution=1.0/16, diameter=1.0 should be at grid index 16
    const diameter_grid_index = @as(usize, @intFromFloat(hs.diameter / hs.resolution));
    try t.expectEqual(@as(usize, 16), diameter_grid_index);
}

test "HardSphereDFT weightFn1" {
    const hs = HardSphereDFT.init(1.0);

    try t.expectApproxEqAbs(0.90724, hs.weightFn1(0), 1e-5);

    // Should be positive deep inside the sphere
    try t.expect(hs.weightFn1(0.5) > 0);

    // Should be negative at contact distance
    try t.expect(hs.weightFn1(1.0) < 0);

    // Test continuity at contact distance
    const eps = 1e-6;
    try t.expectApproxEqAbs(hs.weightFn1(1.0 - eps), hs.weightFn1(1.0 + eps), 1e-3);

    // Test oscillating structure of the tail
    try t.expect(hs.weightFn1(1.5) < 0);
    try t.expect(hs.weightFn1(2.0) > 0);
    try t.expect(hs.weightFn1(2.5) < 0);
}

const WeightIntegral = struct {
    /// The expansion coefficients of density-independent weights over the grid
    /// combinations, or w_i(|r - r'|)
    /// [0]: weightFn0, [1]: weightFn1 (unified), [2]: weightFn2
    kernels: [3]conv.Kernel,

    /// Reference to the grid
    grid: Grid,

    /// Reference to the hard sphere DFT
    hs: HardSphereDFT,

    /// Allocator for memory management
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: Grid, hs: HardSphereDFT) !Self {
        const n = grid.points.len;
        var kernels: [3]conv.Kernel = undefined;
        const weight_fns = hs.getWeightFns();

        var simpson_weights: [3]conv.Weights = undefined;
        for (0..3) |k| {
            const weight_fn = weight_fns[k];
            // Convert relative distances to actual grid indices: distance * GRID_PPD
            var grid_nodes: [3]usize = undefined;
            for (weight_fn.nodes, 0..) |relative_dist, idx| {
                grid_nodes[idx] = relative_dist * GRID_PPD;
            }
            const actual_nodes = grid_nodes[0..weight_fn.nodes.len];

            // Always use Simpson weights since kernel_size is always 2n+1 with GRID_PPD=16
            simpson_weights[k] = try conv.Weights.init(.simpson, allocator, actual_nodes, grid.spacing);
        }
        defer for (0..3) |k| simpson_weights[k].deinit(allocator);

        // Create the three kernels
        for (0..3) |i| {
            const weight_fn = weight_fns[i];
            const max_distance = weight_fn.nodes[weight_fn.nodes.len - 1];
            const kernel_size = max_distance * GRID_PPD + 1;

            var kernel_values = try allocator.alloc(f64, kernel_size);
            defer allocator.free(kernel_values);

            for (0..kernel_size) |offset| {
                const distance = grid.spacing * @as(f64, @floatFromInt(offset));
                kernel_values[offset] = weight_fn.ptr(hs, distance);
            }

            kernels[i] = try conv.Kernel.init(allocator, kernel_values, n, simpson_weights[i]);
        }

        return Self{
            .kernels = kernels,
            .grid = grid,
            .hs = hs,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        // Free weight kernels
        for (0..3) |i| {
            self.kernels[i].deinit(self.allocator);
        }
    }
};

test "WeightIntegral init" {
    const allocator = t.allocator;
    const hs = HardSphereDFT.init(1.0);
    var grid = try Grid.init(allocator, hs.resolution, 5.0);
    defer grid.deinit();
    var kernel = try WeightIntegral.init(allocator, grid, hs);
    defer kernel.deinit();
}

pub const HardSphereWorkspace = struct {
    /// The smoothed (weighted) density profile
    weighted_density: []f64,

    /// The functional derivative of the smoothed density profile with respect to density
    weighted_density_derivative: []f64,

    /// The expansion coefficients of the smoothed (weighted) density profiles (array of 3 vectors)
    weighted_density_expansions: [3][]f64,

    /// Fourier-transformed density profile (complex vector)
    density_fourier: []math.Complex(f64),

    /// Fourier-transformed results for the weighted densities (complex vector)
    weighted_density_fourier: []math.Complex(f64),

    /// Fourier-transformed expansion coefficients
    weight_fns_fourier: [3][]math.Complex(f64),

    /// Allocator used for memory management
    allocator: std.mem.Allocator,

    /// Reference to the grid
    grid: Grid,

    /// Reference to the hard sphere DFT
    hs: HardSphereDFT,

    const Self = @This();

    /// Initialize the density fnal workspace buffers
    pub fn init(allocator: std.mem.Allocator, grid: Grid, hs: HardSphereDFT) !Self {
        var weighted_density_expansions: [3][]f64 = undefined;
        var weight_fns_fourier: [3][]math.Complex(f64) = undefined;

        const n = grid.points.len;
        const fft_size = grid.getFftSize();

        // Allocate memory for each array in weighted_density_expansions
        for (0..3) |i| {
            weighted_density_expansions[i] = try allocator.alloc(f64, n);
        }

        // Allocate memory for Fourier-transformed weights
        for (0..3) |i| {
            weight_fns_fourier[i] = try allocator.alloc(math.Complex(f64), fft_size);
        }

        return Self{
            .weighted_density = try allocator.alloc(f64, n),
            .weighted_density_derivative = try allocator.alloc(f64, n),
            .weighted_density_expansions = weighted_density_expansions,
            .density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weighted_density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weight_fns_fourier = weight_fns_fourier,
            .allocator = allocator,
            .grid = grid,
            .hs = hs,
        };
    }

    /// Deinitialize and free allocated memory
    pub fn deinit(self: *Self) void {
        {
            self.allocator.free(self.weighted_density);
            self.allocator.free(self.weighted_density_derivative);

            for (0..3) |i| {
                self.allocator.free(self.weighted_density_expansions[i]);
            }

            for (0..3) |i| {
                self.allocator.free(self.weight_fns_fourier[i]);
            }

            self.allocator.free(self.density_fourier);
            self.allocator.free(self.weighted_density_fourier);
        }
    }
};

test "HardSphereWorkspace init" {
    const allocator = t.allocator;
    const hs = HardSphereDFT.init(1.0);
    var grid = try Grid.init(allocator, hs.resolution, 5.0);
    defer grid.deinit();
    var workspace = try HardSphereWorkspace.init(allocator, grid, hs);
    defer workspace.deinit();
}
