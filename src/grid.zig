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

    /// Initialize the grid with specified spacing, expanding to fit min_size
    pub fn init(allocator: std.mem.Allocator, spacing: f64, min_size: f64) !Self {
        // Calculate number of points needed to fit at least min_size with given spacing
        const n = @as(usize, @intFromFloat(@ceil(min_size / spacing))) + 1;
        const actual_size = @as(f64, @floatFromInt(n - 1)) * spacing;

        const points = try allocator.alloc(f64, n);

        // Pre-compute all grid points
        for (0..n) |i| {
            points[i] = @as(f64, @floatFromInt(i)) * spacing;
        }

        return Self{
            .spacing = spacing,
            .size = actual_size,
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

    /// Get distance between two grid points by their indices
    pub fn distance(self: Self, i: usize, j: usize) f64 {
        return @abs(self.points[i] - self.points[j]);
    }
};

test "Grid init" {
    var grid = try Grid.init(t.allocator, 0.1, 10.0);
    defer grid.deinit();

    try t.expectEqual(@as(usize, 101), grid.points.len);
    try t.expectEqual(@as(f64, 10.0), grid.size);
    try t.expectEqual(@as(f64, 0.1), grid.spacing);
}

test "Grid points array" {
    var grid = try Grid.init(t.allocator, 0.5, 5.0);
    defer grid.deinit();

    try t.expectEqual(@as(f64, 0.0), grid.points[0]);
    try t.expectEqual(@as(f64, 0.5), grid.points[1]);
    try t.expectEqual(@as(f64, 5.0), grid.points[10]);
}

test "Grid FFT size" {
    var grid = try Grid.init(t.allocator, 0.5, 5.0);
    defer grid.deinit();

    try t.expectEqual(@as(usize, 32), grid.getFftSize());
}

test "Grid distance" {
    var grid = try Grid.init(t.allocator, 0.5, 5.0);
    defer grid.deinit();

    try t.expectEqual(@as(f64, 0.0), grid.distance(0, 0));
    try t.expectEqual(@as(f64, 0.5), grid.distance(0, 1));
    try t.expectEqual(@as(f64, 0.5), grid.distance(1, 0));
    try t.expectEqual(@as(f64, 5.0), grid.distance(0, 10));
}

test "Grid expansion with significant remainder" {
    // Test case where min_size / spacing = 5.0 / 0.3 = 16.667...
    // Should expand to accommodate the spacing
    var grid = try Grid.init(t.allocator, 0.3, 5.0);
    defer grid.deinit();

    // Expected: n = ceil(5.0/0.3) + 1 = ceil(16.667) + 1 = 17 + 1 = 18 points
    // actual_size = (18-1) * 0.3 = 17 * 0.3 = 5.1
    try t.expectEqual(@as(usize, 18), grid.points.len);
    try t.expectEqual(@as(f64, 5.1), grid.size);
    try t.expectEqual(@as(f64, 0.3), grid.spacing);

    // Verify the grid actually covers the min_size and more
    try t.expect(grid.size >= 5.0);

    // Check specific grid points
    try t.expectEqual(@as(f64, 0.0), grid.points[0]);
    try t.expectEqual(@as(f64, 0.3), grid.points[1]);
    try t.expectEqual(@as(f64, 5.1), grid.points[17]); // Last point
}
