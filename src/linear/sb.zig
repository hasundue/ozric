const std = @import("std");
const testing = std.testing;
const Allocator = std.mem.Allocator;

/// Symmetric band matrix structure
pub const SymmetricBandMatrix = struct {
    data: []f64,
    n: usize,
    k: usize,
    lda: usize,

    const Self = @This();

    /// Initialize a symmetric band matrix with given dimensions
    pub fn init(allocator: Allocator, n: usize, k: usize) !Self {
        const lda = k + 1;
        const data = try allocator.alloc(f64, lda * n);

        return Self{
            .data = data,
            .n = n,
            .k = k,
            .lda = lda,
        };
    }

    /// Deinitialize and free the matrix data
    pub fn deinit(self: Self, allocator: Allocator) void {
        allocator.free(self.data);
    }

    /// Clear the matrix (set all elements to zero)
    pub fn clear(self: *Self) void {
        @memset(self.data, 0.0);
    }

    /// Set matrix element at (i, j) using upper triangle storage
    /// Only sets elements within the band (|i - j| <= k)
    pub fn set(self: *Self, i: usize, j: usize, value: f64) void {
        if (i >= self.n or j >= self.n) return;

        // Only store upper triangle elements (i <= j)
        if (i > j) return;

        // Check if element is within the band
        if (j - i > self.k) return;

        const band_row = self.k + i - j;
        const idx = band_row + j * self.lda;

        if (idx < self.data.len) {
            self.data[idx] = value;
        }
    }

    /// Get matrix element at (i, j) using symmetry
    pub fn get(self: Self, i: usize, j: usize) f64 {
        if (i >= self.n or j >= self.n) return 0.0;

        // Use symmetry: A(i,j) = A(j,i)
        const row = @min(i, j);
        const col = @max(i, j);

        // Check if element is within the band
        if (col - row > self.k) return 0.0;

        const band_row = self.k + row - col;
        const idx = band_row + col * self.lda;

        if (idx < self.data.len) {
            return self.data[idx];
        }

        return 0.0;
    }
};

test "symmetric band matrix basic operations" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    // Create a 3x3 matrix with bandwidth 1
    var matrix = try SymmetricBandMatrix.init(allocator, 3, 1);
    defer matrix.deinit(allocator);

    // Set some values
    matrix.set(0, 0, 2.0);
    matrix.set(0, 1, 1.0);
    matrix.set(1, 1, 2.0);
    matrix.set(1, 2, 1.0);
    matrix.set(2, 2, 2.0);

    // Test getting values with symmetry
    const diagonal_positions = [_][2]usize{ .{ 0, 0 }, .{ 1, 1 }, .{ 2, 2 } };
    const diagonal_expected = [_]f64{ 2.0, 2.0, 2.0 };
    for (diagonal_positions, diagonal_expected) |pos, expected| {
        try eqa(expected, matrix.get(pos[0], pos[1]), 1e-10);
    }

    // Test off-diagonal elements and symmetry
    const off_diag_positions = [_][2]usize{ .{ 0, 1 }, .{ 1, 0 }, .{ 1, 2 }, .{ 2, 1 } };
    const off_diag_expected = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    for (off_diag_positions, off_diag_expected) |pos, expected| {
        try eqa(expected, matrix.get(pos[0], pos[1]), 1e-10);
    }

    // Test out-of-band elements
    const out_of_band_positions = [_][2]usize{ .{ 0, 2 }, .{ 2, 0 } };
    for (out_of_band_positions) |pos| {
        try eqa(0.0, matrix.get(pos[0], pos[1]), 1e-10);
    }
}

test "symmetric band matrix properties" {
    const allocator = testing.allocator;
    const eq = testing.expectEqual;

    var matrix = try SymmetricBandMatrix.init(allocator, 5, 2);
    defer matrix.deinit(allocator);

    try eq(5, matrix.n);
    try eq(2, matrix.k);
    try eq(3, matrix.lda);
    try eq(15, matrix.data.len);
}

test "matrix clear operation" {
    const allocator = testing.allocator;

    var matrix = try SymmetricBandMatrix.init(allocator, 3, 1);
    defer matrix.deinit(allocator);

    // Set some values
    matrix.set(0, 0, 5.0);
    matrix.set(0, 1, 3.0);

    // Clear the matrix
    matrix.clear();

    // Check that entire data array is zeroed
    for (matrix.data) |val| {
        try testing.expectEqual(@as(f64, 0.0), val);
    }
}
