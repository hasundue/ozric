const std = @import("std");
const c = @cImport({
    @cInclude("cblas.h");
});

/// Grid parameters for radial functions
pub const GridParams = struct {
    n_points: usize,
    r_max: f64,
    dr: f64,

    pub fn init(n_points: usize, r_max: f64) GridParams {
        return GridParams{
            .n_points = n_points,
            .r_max = r_max,
            .dr = r_max / @as(f64, @floatFromInt(n_points - 1)),
        };
    }
};

/// Radial function data structure
pub const RadialFunction = struct {
    values: []f64,
    grid: GridParams,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: GridParams) !Self {
        const values = try allocator.alloc(f64, grid.n_points);
        @memset(values, 0.0);
        return Self{
            .values = values,
            .grid = grid,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.values);
    }

    pub fn getRadius(self: *const Self, i: usize) f64 {
        return @as(f64, @floatFromInt(i)) * self.grid.dr;
    }
};

/// Vector for linear algebra operations with BLAS
pub const Vector = struct {
    data: []f64,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, size: usize) !Self {
        const data = try allocator.alloc(f64, size);
        return Self{
            .data = data,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.data);
    }

    /// Compute dot product using BLAS
    pub fn dot(self: *const Self, other: *const Self) f64 {
        std.debug.assert(self.data.len == other.data.len);
        const n: c_int = @intCast(self.data.len);
        return c.cblas_ddot(n, self.data.ptr, 1, other.data.ptr, 1);
    }

    /// Scale vector by a constant using BLAS
    pub fn scale(self: *Self, alpha: f64) void {
        const n: c_int = @intCast(self.data.len);
        c.cblas_dscal(n, alpha, self.data.ptr, 1);
    }
};
