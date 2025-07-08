const std = @import("std");
const math = std.math;
const t = std.testing;
const util = @import("util.zig");

pub const HardSphereDensityFunctional = struct {
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

    /// The expansion coefficients of density-independent weights over the grid
    /// combinations, or weight_functions(|r - r'|) - 3 matrices of size N×N
    weight_functions: [3][][]f64,

    /// Fourier-transformed expansion coefficients
    weight_functions_fourier: [3][]math.Complex(f64),

    /// Grid points (radial positions)
    grid_points: []f64,

    /// Grid spacing (radial step size)
    grid_spacing: f64,

    /// System size
    system_size: f64,

    /// Number of grid points
    num_grid_points: usize,

    /// Flag to check if buffers are initialized
    initialized: bool,

    /// Flag to check if weights are computed
    weights_computed: bool,

    /// Allocator used for memory management
    allocator: std.mem.Allocator,

    const Self = @This();

    /// Initialize the density functional buffers
    pub fn init(allocator: std.mem.Allocator, num_grid_points: usize, system_size: f64) !Self {
        const fft_size = util.nextPowerOfTwo(2 * num_grid_points);
        var weighted_density_expansions: [3][]f64 = undefined;
        var weight_functions: [3][][]f64 = undefined;
        var weight_functions_fourier: [3][]math.Complex(f64) = undefined;

        // Allocate memory for each array in weighted_density_expansions
        for (0..3) |i| {
            weighted_density_expansions[i] = try allocator.alloc(f64, num_grid_points);
        }

        // Allocate memory for weight matrices (N×N)
        for (0..3) |i| {
            weight_functions[i] = try allocator.alloc([]f64, num_grid_points);
            for (0..num_grid_points) |j| {
                weight_functions[i][j] = try allocator.alloc(f64, num_grid_points);
            }
        }

        // Allocate memory for Fourier-transformed weights
        for (0..3) |i| {
            weight_functions_fourier[i] = try allocator.alloc(math.Complex(f64), fft_size);
        }

        return Self{
            .weighted_density = try allocator.alloc(f64, num_grid_points),
            .weighted_density_derivative = try allocator.alloc(f64, num_grid_points),
            .weighted_density_expansions = weighted_density_expansions,
            .density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weighted_density_fourier = try allocator.alloc(math.Complex(f64), fft_size),
            .weight_functions = weight_functions,
            .weight_functions_fourier = weight_functions_fourier,
            .grid_points = try allocator.alloc(f64, num_grid_points),
            .grid_spacing = system_size / @as(f64, @floatFromInt(num_grid_points)),
            .system_size = system_size,
            .num_grid_points = num_grid_points,
            .initialized = true,
            .weights_computed = false,
            .allocator = allocator,
        };
    }

    /// Deinitialize and free allocated memory
    pub fn deinit(self: *Self) void {
        if (self.initialized) {
            self.allocator.free(self.weighted_density);
            self.allocator.free(self.weighted_density_derivative);
            self.allocator.free(self.grid_points);

            for (0..3) |i| {
                self.allocator.free(self.weighted_density_expansions[i]);
            }

            // Free weight matrices
            for (0..3) |i| {
                for (0..self.num_grid_points) |j| {
                    self.allocator.free(self.weight_functions[i][j]);
                }
                self.allocator.free(self.weight_functions[i]);
            }

            // Free Fourier-transformed weights
            for (0..3) |i| {
                self.allocator.free(self.weight_functions_fourier[i]);
            }

            self.allocator.free(self.density_fourier);
            self.allocator.free(self.weighted_density_fourier);

            self.initialized = false;
            self.weights_computed = false;
        }
    }

    /// Initialize the grid points
    pub fn initializeGrid(self: *Self) void {
        for (0..self.num_grid_points) |i| {
            self.grid_points[i] = @as(f64, @floatFromInt(i)) * self.grid_spacing;
        }
    }

    /// Precompute w_i_ over the grid combinations
    pub fn precomputeWeights(self: *Self, sigma: f64) void {
        // Skip calculation if weights are already computed
        if (self.weights_computed) {
            return;
        }

        // Compute weights for upper triangle (i <= j)
        for (0..self.num_grid_points) |i| {
            for (i..self.num_grid_points) |j| {
                const r_distance = @abs(self.grid_points[i] - self.grid_points[j]);
                self.weight_functions[0][i][j] = weightFunctionZeroth(sigma, r_distance);
                self.weight_functions[1][i][j] = weightFunctionFirst(sigma, r_distance);
                self.weight_functions[2][i][j] = weightFunctionSecond(sigma, r_distance);
            }
        }

        // Fill lower triangle by symmetry
        for (0..self.num_grid_points) |i| {
            for (0..i) |j| {
                self.weight_functions[0][i][j] = self.weight_functions[0][j][i];
                self.weight_functions[1][i][j] = self.weight_functions[1][j][i];
                self.weight_functions[2][i][j] = self.weight_functions[2][j][i];
            }
        }

        // Mark weights as computed
        self.weights_computed = true;
    }
};

pub fn weightFunctionZeroth(
    sigma: f64,
    r: f64,
) f64 {
    if (r > sigma) return 0.0;

    const sigma3 = math.pow(f64, sigma, 3);
    const pi = math.pi;

    return 3 / (4 * pi * sigma3);
}

pub fn weightFunctionFirst(
    sigma: f64,
    r: f64,
) f64 {
    const x = r / sigma;
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

test weightFunctionFirst {
    const sigma = 1.0;

    const w1 = struct {
        pub fn call(r: f64) f64 {
            return weightFunctionFirst(sigma, r);
        }
    }.call;

    try t.expectApproxEqAbs(0.90724, w1(0), 1e-5);

    // Should be positive deep inside the sphere
    try t.expect(w1(0.5) > 0);

    // Should be negative at contact distance
    try t.expect(w1(1.0) < 0);

    // Test continuity at contact distance
    const eps = 1e-6;
    try t.expectApproxEqAbs(w1(sigma - eps), w1(sigma + eps), 1e-3);

    // Test oscillating structure of the tail
    try t.expect(w1(1.5) < 0);
    try t.expect(w1(2.0) > 0);
    try t.expect(w1(2.5) < 0);
}

pub fn weightFunctionSecond(
    sigma: f64,
    r: f64,
) f64 {
    if (r >= sigma) return 0.0;

    const x = r / sigma;
    const x2 = x * x;
    const pi = math.pi;

    return (15 / 8 / pi / sigma) * (1 - 3 * x + 3 * x2);
}
