const std = @import("std");
const math = std.math;
const t = std.testing;
const conv = @import("convolution.zig");

const Grid = @import("grid.zig").Grid;
const Matrix = @import("matrix.zig").Matrix;

pub const HardSphereDFT = struct {
    /// Hard-sphere diameter
    diameter: f64,

    const Self = @This();

    pub fn init(diameter: f64) Self {
        return Self{
            .diameter = diameter,
        };
    }

    pub fn weightFn0(self: Self, r: f64) f64 {
        if (r > self.diameter) return 0.0;

        const sigma3 = math.pow(f64, self.diameter, 3);
        const pi = math.pi;

        return 3 / (4 * pi * sigma3);
    }

    pub fn weightFn1(self: Self, r: f64) f64 {
        if (r < self.diameter) {
            return self.weightFn1c(r);
        } else {
            return self.weightFn1t(r);
        }
    }

    // The first order weight function inside the hard sphere
    pub fn weightFn1c(self: Self, r: f64) f64 {
        if (r > self.diameter) return 0.0;

        const x = r / self.diameter;
        const x2 = math.pow(f64, x, 2);

        const a_0 = 0.90724;
        const a_1 = -1.23717;
        const a_2 = 0.21616;

        return a_0 + a_1 * x + a_2 * x2;
    }

    /// The second order weight function outside the hard sphere
    pub fn weightFn1t(self: Self, r: f64) f64 {
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

pub const HardSphereKernels = struct {
    /// The expansion coefficients of density-independent weights over the grid
    /// combinations, or w_i(|r - r'|)
    weights: [4]conv.Kernel,

    /// Reference to the grid
    grid: Grid,

    /// Reference to the hard sphere DFT
    hs: HardSphereDFT,

    /// Allocator for memory management
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: Grid, hs: HardSphereDFT) !Self {
        const n = grid.points.len;
        var weight_fns: [3]Matrix = undefined;

        for (0..3) |i| {
            weight_fns[i] = try Matrix.init(allocator, n, n);
        }

        // Fill in the upper triangle of the matrices
        for (0..n) |i| {
            for (i..n) |j| {
                const d = grid.distance(i, j);
                weight_fns[0].ptr(i, j).* = hs.weightFn0(d);
                weight_fns[1].ptr(i, j).* = hs.weightFn1(d);
                weight_fns[2].ptr(i, j).* = hs.weightFn2(d);
            }
        }

        // Fill in the lower triangle of the matrices by symmetry
        for (0..n) |i| {
            for (0..i) |j| {
                for (0..3) |k| {
                    weight_fns[k].ptr(i, j).* = weight_fns[k].at(j, i);
                }
            }
        }

        return Self{
            .weights = weight_fns,
            .grid = grid,
            .hs = hs,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        // Free weight matrices
        for (0..3) |i| {
            self.weights[i].deinit(self.allocator);
        }
    }
};

test "HardSphereKernel init" {
    const allocator = t.allocator;
    var grid = try Grid.init(allocator, 10, 5.0);
    defer grid.deinit();
    const hs = HardSphereDFT.init(1.0);
    var kernel = try HardSphereKernels.init(allocator, grid, hs);
    defer kernel.deinit();

    // Check that the weight functions are initialized correctly
    try t.expect(kernel.weights[0].at(0, 1) == hs.weightFn0(grid.distance(0, 1)));
    try t.expect(kernel.weights[1].at(0, 1) == hs.weightFn1(grid.distance(0, 1)));
    try t.expect(kernel.weights[2].at(0, 1) == hs.weightFn2(grid.distance(0, 1)));
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
    var grid = try Grid.init(allocator, 10, 5.0);
    defer grid.deinit();
    const hs = HardSphereDFT.init(1.0);
    var workspace = try HardSphereWorkspace.init(allocator, grid, hs);
    defer workspace.deinit();
}
