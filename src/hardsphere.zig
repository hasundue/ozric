const std = @import("std");
const math = std.math;
const t = std.testing;
const grid_mod = @import("grid.zig");
const Grid = grid_mod.Grid;

pub const HardSphereDFT = struct {
    /// Hard sphere diameter (sigma)
    diameter: f64,

    const Self = @This();

    /// Initialize the hard sphere system
    pub fn init(diameter: f64) Self {
        return Self{
            .diameter = diameter,
        };
    }

    /// Zeroth-order weight function
    pub fn weightFunctionZeroth(self: Self, r: f64) f64 {
        if (r > self.diameter) return 0.0;

        const sigma3 = math.pow(f64, self.diameter, 3);
        const pi = math.pi;

        return 3 / (4 * pi * sigma3);
    }

    /// First-order weight function
    pub fn weightFunctionFirst(self: Self, r: f64) f64 {
        const x = r / self.diameter;
        const x2 = math.pow(f64, x, 2);

        if (x <= 1.0) {
            const a_0 = 0.90724;
            const a_1 = -1.23717;
            const a_2 = 0.21616;

            return a_0 + a_1 * x + a_2 * x2;
        } else {
            const x3 = math.pow(f64, x, 3);

            const c = -0.10244;
            const b_0 = 35.134;
            const alpha = 4.934;
            const b_1 = -98.684;
            const beta_1 = 3.5621;
            const b_2 = 92.693;
            const beta_2 = 12.0;
            const b_3 = -29.257;

            const term1 = c * @exp(-beta_1 * (x - 1.0)) * @sin(alpha * (x - 1.0));
            const term2 = @exp(-beta_2 * (x - 1.0)) * (b_0 + b_1 * x + b_2 * x2 + b_3 * x3);
            return term1 + term2;
        }
    }

    /// Second-order weight function
    pub fn weightFunctionSecond(self: Self, r: f64) f64 {
        if (r >= self.diameter) return 0.0;

        const x = r / self.diameter;
        const x2 = x * x;
        const pi = math.pi;

        return (15 / 8 / pi / self.diameter) * (1 - 3 * x + 3 * x2);
    }
};

test "HardSphereDFT init" {
    const hs = HardSphereDFT.init(1.0);
    try t.expectEqual(@as(f64, 1.0), hs.diameter);
}

test "HardSphereDFT weightFunctionFirst" {
    const hs = HardSphereDFT.init(1.0);

    try t.expectApproxEqAbs(0.90724, hs.weightFunctionFirst(0), 1e-5);

    // Should be positive deep inside the sphere
    try t.expect(hs.weightFunctionFirst(0.5) > 0);

    // Should be negative at contact distance
    try t.expect(hs.weightFunctionFirst(1.0) < 0);

    // Test continuity at contact distance
    const eps = 1e-6;
    try t.expectApproxEqAbs(hs.weightFunctionFirst(1.0 - eps), hs.weightFunctionFirst(1.0 + eps), 1e-3);

    // Test oscillating structure of the tail
    try t.expect(hs.weightFunctionFirst(1.5) < 0);
    try t.expect(hs.weightFunctionFirst(2.0) > 0);
    try t.expect(hs.weightFunctionFirst(2.5) < 0);
}

pub const HardSphereKernel = struct {
    /// The expansion coefficients of density-independent weights over the grid
    /// combinations, or weight_functions(|r - r'|) - 3 matrices of size N×N
    weight_functions: [3][][]f64,

    /// Reference to the grid
    grid: *const Grid,

    /// Reference to the hard sphere DFT
    hs: *const HardSphereDFT,

    /// Allocator for memory management
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: *const Grid, hs: *const HardSphereDFT) !Self {
        const n = grid.points.len;
        var weight_functions: [3][][]f64 = undefined;

        // Allocate memory for weight matrices (N×N)
        for (0..3) |i| {
            weight_functions[i] = try allocator.alloc([]f64, n);
            for (0..n) |j| {
                weight_functions[i][j] = try allocator.alloc(f64, n);
            }
        }

        return Self{
            .weight_functions = weight_functions,
            .grid = grid,
            .hs = hs,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        // Free weight matrices
        for (0..3) |i| {
            for (0..self.grid.points.len) |j| {
                self.allocator.free(self.weight_functions[i][j]);
            }
            self.allocator.free(self.weight_functions[i]);
        }
    }

    pub fn computeWeights(self: *Self) void {
        const n = self.grid.points.len;
        // Compute weights for upper triangle (i <= j)
        for (0..n) |i| {
            for (i..n) |j| {
                const r_distance = @abs(self.grid.points[i] - self.grid.points[j]);
                self.weight_functions[0][i][j] = self.hs.weightFunctionZeroth(r_distance);
                self.weight_functions[1][i][j] = self.hs.weightFunctionFirst(r_distance);
                self.weight_functions[2][i][j] = self.hs.weightFunctionSecond(r_distance);
            }
        }
        // Fill lower triangle by symmetry
        for (0..n) |j| {
            for (j..n) |i| {
                self.weight_functions[0][i][j] = self.weight_functions[0][j][i];
                self.weight_functions[1][i][j] = self.weight_functions[1][j][i];
                self.weight_functions[2][i][j] = self.weight_functions[2][j][i];
            }
        }
    }
};

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
    weight_functions_fourier: [3][]math.Complex(f64),

    /// Allocator used for memory management
    allocator: std.mem.Allocator,

    /// Reference to the grid
    grid: *const Grid,

    /// Reference to the hard sphere DFT
    hs: *const HardSphereDFT,

    const Self = @This();

    /// Initialize the density functional workspace buffers
    pub fn init(allocator: std.mem.Allocator, grid: *const Grid, hs: *const HardSphereDFT) !Self {
        var weighted_density_expansions: [3][]f64 = undefined;
        var weight_functions_fourier: [3][]math.Complex(f64) = undefined;

        const n = grid.points.len;
        const fft_size = grid.getFftSize();

        // Allocate memory for each array in weighted_density_expansions
        for (0..3) |i| {
            weighted_density_expansions[i] = try allocator.alloc(f64, n);
        }

        // Allocate memory for Fourier-transformed weights
        for (0..3) |i| {
            weight_functions_fourier[i] = try allocator.alloc(math.Complex(f64), fft_size);
        }

        return Self{
            .weighted_density = try allocator.alloc(f64, n),
            .weighted_density_derivative = try allocator.alloc(f64, n),
            .weighted_density_expansions = weighted_density_expansions,
            .density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weighted_density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weight_functions_fourier = weight_functions_fourier,
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

            // Free Fourier-transformed weights
            for (0..3) |i| {
                self.allocator.free(self.weight_functions_fourier[i]);
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
    var workspace = try HardSphereWorkspace.init(allocator, &grid, &hs);
    defer workspace.deinit();
}
