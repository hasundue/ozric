const std = @import("std");

pub const Matrix = struct {
    data: []f64,
    rows: usize,
    cols: usize,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) !Matrix {
        const data = try allocator.alloc(f64, rows * cols);
        return Matrix{
            .data = data,
            .rows = rows,
            .cols = cols,
        };
    }

    pub fn deinit(self: *Self, allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }

    pub fn at(self: Self, row: usize, col: usize) f64 {
        return self.data[row * self.cols + col];
    }

    pub fn ptr(self: *Self, row: usize, col: usize) *f64 {
        return &self.data[row * self.cols + col];
    }

    // pub fn row(self: Self, row: usize) []f64 {
    //     return self.data[row * self.cols ..][0..self.cols];
    // }
};
