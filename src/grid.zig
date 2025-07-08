const std = @import("std");
const util = @import("util.zig");
const t = std.testing;

pub const Grid = struct {
    /// Grid spacing (radial step size)
    spacing: f64,

    /// System size
    size: f64,

    /// Pre-computed grid points array
    points: []f64,

    /// Allocator for memory management
    allocator: std.mem.Allocator,

    const Self = @This();

    /// Initialize the grid with pre-computed points
    pub fn init(allocator: std.mem.Allocator, n: usize, size: f64) !Self {
        const spacing = size / @as(f64, @floatFromInt(n));
        const points = try allocator.alloc(f64, n);

        // Pre-compute all grid points
        for (0..n) |i| {
            points[i] = @as(f64, @floatFromInt(i)) * spacing;
        }

        return Self{
            .spacing = spacing,
            .size = size,
            .points = points,
            .allocator = allocator,
        };
    }

    /// Clean up allocated memory
    pub fn deinit(self: *Self) void {
        self.allocator.free(self.points);
    }

    /// Get FFT size for calculations
    pub fn getFftSize(self: Self) usize {
        return util.nextPowerOfTwo(2 * self.points.len);
    }
};

test "Grid init" {
    var grid = try Grid.init(t.allocator, 100, 10.0);
    defer grid.deinit();

    try t.expectEqual(@as(usize, 100), grid.points.len);
    try t.expectEqual(@as(f64, 10.0), grid.size);
    try t.expectEqual(@as(f64, 0.1), grid.spacing);
    try t.expectEqual(@as(usize, 100), grid.points.len);
}

test "Grid points array" {
    var grid = try Grid.init(t.allocator, 10, 5.0);
    defer grid.deinit();

    try t.expectEqual(@as(f64, 0.0), grid.points[0]);
    try t.expectEqual(@as(f64, 0.5), grid.points[1]);
    try t.expectEqual(@as(f64, 4.5), grid.points[9]);
}

test "Grid FFT size" {
    var grid = try Grid.init(t.allocator, 10, 5.0);
    defer grid.deinit();

    try t.expectEqual(@as(usize, 32), grid.getFftSize());
}
