const std = @import("std");
const math = std.math;
const t = std.testing;
const conv = @import("convolution.zig");

const Grid = @import("grid.zig").Grid;

pub const HardSphereOptions = struct {
    min_resolution: f64,
};

pub const HardSphereDFT = struct {
    /// Hard-sphere diameter
    diameter: f64,

    /// Pre-calculated optimal grid spacing that aligns diameter with a grid point
    resolution: f64,

    /// Relative diameter of hard sphere to the resolution
    size: usize,

    /// Descritized weight functions over the grid
    weight_functions: [3][]f64,

    const Self = @This();

    /// Number of grid points per diameter for optimal discretization
    const GRID_PPD = 16;

    pub fn init(
        allocator: std.mem.Allocator,
        diameter: f64,
    ) !Self {
        const size = GRID_PPD;
        const resolution: f64 = diameter / @as(f64, size);

        const weights = .{ weightFn0, weightFn1, weightFn2 };
        var weight_functions: [3][]f64 = undefined;
        inline for (0..3) |i| {
            const support = if (i == 1) 2 else 1;
            weight_functions[i] = try allocator.alloc(f64, size * support + 1);
            for (0..weight_functions[i].len) |j| {
                const r = @as(f64, @floatFromInt(j)) * resolution;
                weight_functions[i][j] = weights[i](diameter, r);
            }
        }

        return Self{
            .diameter = diameter,
            .resolution = resolution,
            .size = size,
            .weight_functions = weight_functions,
        };
    }

    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        for (0..3) |i| {
            allocator.free(self.weight_functions[i]);
        }
    }
};

test "HardSphereDFT init" {
    const hs = try HardSphereDFT.init(t.allocator, 1.0);
    defer hs.deinit(t.allocator);

    try t.expectEqual(1.0, hs.diameter);
    try t.expectEqual(1.0 / 16.0, hs.resolution);
    try t.expectEqual(16, hs.size);
    try t.expectEqual(16 + 1, hs.weight_functions[0].len);
    try t.expectEqual(32 + 1, hs.weight_functions[1].len);
    try t.expectEqual(16 + 1, hs.weight_functions[2].len);
}

const WeightIntegral = struct {
    /// The expansion coefficients of density-independent weights over the grid
    /// combinations, or w_i(|r - r'|), weighted by integration weights
    kernels: [3]conv.Kernel,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, hs: HardSphereDFT, grid: Grid) !Self {
        var kernels: [3]conv.Kernel = undefined;
        for (0..3) |i| {
            const simpson_weights = try conv.RadialWeights.init(
                .simpson,
                allocator,
                hs.weight_functions[i].len,
                grid.spacing,
            );
            defer simpson_weights.deinit(allocator);

            kernels[i] = try conv.Kernel.init(
                allocator,
                hs.weight_functions[i],
                grid.points.len,
                simpson_weights,
            );
        }
        return Self{ .kernels = kernels };
    }

    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        // Free weight kernels
        for (0..3) |i| {
            self.kernels[i].deinit(allocator);
        }
    }
};

test "WeightIntegral init" {
    const allocator = t.allocator;

    const hs = try HardSphereDFT.init(allocator, 1.0);
    defer hs.deinit(allocator);

    var grid = try Grid.init(allocator, hs.resolution, 5.0);
    defer grid.deinit();

    var kernel = try WeightIntegral.init(allocator, hs, grid);
    defer kernel.deinit(allocator);
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

    const hs = try HardSphereDFT.init(allocator, 1.0);
    defer hs.deinit(allocator);

    var grid = try Grid.init(allocator, hs.resolution, 5.0);
    defer grid.deinit();

    var workspace = try HardSphereWorkspace.init(allocator, grid, hs);
    defer workspace.deinit();
}

fn weightFn0(sigma: f64, r: f64) f64 {
    if (r > sigma) return 0.0;

    const sigma3 = math.pow(f64, sigma, 3);
    const pi = math.pi;

    return 3 / (4 * pi * sigma3);
}

fn weightFn1(sigma: f64, r: f64) f64 {
    const x = r / sigma;
    const x2 = math.pow(f64, x, 2);
    const x3 = math.pow(f64, x, 3);

    if (r <= sigma) {
        const a_0 = 0.90724;
        const a_1 = -1.23717;
        const a_2 = 0.21616;

        return a_0 + a_1 * x + a_2 * x2;
    } else {
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
}

test "weightFn1" {
    try t.expectApproxEqAbs(0.90724, weightFn1(1.0, 0), 1e-5);

    // Should be positive deep inside the sphere
    try t.expect(weightFn1(1.0, 0.5) > 0);

    // Should be negative at contact distance
    try t.expect(weightFn1(1.0, 1.0) < 0);

    // Test continuity at contact distance
    const eps = 1e-6;
    try t.expectApproxEqAbs(weightFn1(1.0, 1.0 - eps), weightFn1(1.0, 1.0 + eps), 1e-3);

    // Test oscillating structure of the tail
    try t.expect(weightFn1(1.0, 1.5) < 0);
    try t.expect(weightFn1(1.0, 2.0) > 0);
    try t.expect(weightFn1(1.0, 2.5) < 0);
}

fn weightFn2(sigma: f64, r: f64) f64 {
    if (r > sigma) return 0.0;

    const x = r / sigma;
    const x2 = x * x;
    const pi: f64 = math.pi;

    return 5 / (4 * pi) * (6 - 12 * x + 5 * x2);
}
