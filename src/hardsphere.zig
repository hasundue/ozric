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
    /// [0] = w_3(z), [1] = w_2(z), [2] = w^i_2(z)
    weights: [3][]f64,

    const Self = @This();

    /// Number of grid points per diameter for optimal discretization
    const GRID_PPD = 16;

    pub fn init(
        allocator: std.mem.Allocator,
        diameter: f64,
    ) !Self {
        const size = GRID_PPD;
        const resolution: f64 = diameter / @as(f64, size);

        var weights: [3][]f64 = undefined;
        inline for (0..3) |i| {
            weights[i] = try allocator.alloc(f64, size + 1);
        }
        for (0..size + 1) |i| {
            const i_ = @as(f64, @floatFromInt(i));
            const z = i_ * resolution;
            const x = i_ / @as(f64, size);
            weights[0][i] = math.pi * diameter * diameter * (1.0 - x * x);
            weights[2][i] = 2.0 * math.pi * z;
        }
        @memset(weights[1], 2.0 * math.pi * diameter);

        return Self{
            .diameter = diameter,
            .resolution = resolution,
            .size = size,
            .weights = weights,
        };
    }

    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        for (0..3) |i| {
            allocator.free(self.weights[i]);
        }
    }
};

test "HardSphereDFT.init" {
    const hs = try HardSphereDFT.init(t.allocator, 1.0);
    defer hs.deinit(t.allocator);

    try t.expectEqual(1.0, hs.diameter);
    try t.expectEqual(1.0 / 16.0, hs.resolution);
    try t.expectEqual(16, hs.size);
    try t.expectEqual(16 + 1, hs.weights[0].len);
    try t.expectEqual(16 + 1, hs.weights[1].len);
    try t.expectEqual(16 + 1, hs.weights[2].len);
}

const Kernels = struct {
    /// Convolution kernels for the weight functions
    /// [0] = w3(z), [1] = w2(z), [2] = v2(z)
    weights: [3]conv.RadialKernel,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, hs: HardSphereDFT, grid: Grid) !Self {
        var weights: [3]conv.RadialKernel = undefined;
        for (0..3) |i| {
            const simpson_weights = try conv.RadialWeights.init(
                .simpson,
                allocator,
                hs.weights[i].len,
                grid.spacing,
            );
            defer simpson_weights.deinit(allocator);

            weights[i] = try conv.RadialKernel.init(
                allocator,
                hs.weights[i],
                grid.points.len,
                simpson_weights,
            );
        }
        return Self{ .weights = weights };
    }

    pub fn deinit(self: Self, allocator: std.mem.Allocator) void {
        for (0..3) |i| {
            self.weights[i].deinit(allocator);
        }
    }
};

test "Kernels.init" {
    const allocator = t.allocator;

    const hs = try HardSphereDFT.init(allocator, 1.0);
    defer hs.deinit(allocator);

    var grid = try Grid.init(allocator, hs.resolution, 5.0);
    defer grid.deinit();

    var kernel = try Kernels.init(allocator, hs, grid);
    defer kernel.deinit(allocator);
}

const WeightedDensity = struct {
    data: []f64,
    grid: Grid,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: Grid) !Self {
        const n = grid.points.len;
        return Self{
            .data = try allocator.alloc(f64, 3 * n),
            .grid = grid,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }

    pub fn update(
        self: Self,
        density: []const f64,
        kernels: Kernels,
    ) void {
        const n = density.len;
        kernels.weights[0].convolve(density, self.data[0..n]);
        kernels.weights[1].convolve(density, self.data[n .. 2 * n]);
        kernels.weights[2].convolve(density, self.data[2 * n .. 3 * n]);
    }
};

test "WeightedDensity" {
    const allocator = t.allocator;

    const dft = try HardSphereDFT.init(allocator, 1.0);
    defer dft.deinit(allocator);

    const grid = try Grid.init(allocator, dft.resolution, 5.0);
    defer grid.deinit();
    const n = grid.points.len;
    const center = grid.points.len / 2 - 1;

    const integral = try Kernels.init(allocator, dft, grid);
    defer integral.deinit(allocator);

    var weighted = try WeightedDensity.init(allocator, grid);
    defer weighted.deinit(allocator);

    const density = try allocator.alloc(f64, grid.points.len);
    defer allocator.free(density);
    @memset(density, 1.0);

    weighted.update(density, integral);

    std.debug.print("w[0]: {}\n", .{weighted.data[center]});
    std.debug.print("w[1]: {}\n", .{weighted.data[n + center]});
    std.debug.print("w[2]: {}\n", .{weighted.data[2 * n + center]});
}
