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
    pub fn deinit(self: *Self, allocator: Allocator) void {
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

/// Create a symmetric band matrix from a kernel for convolution
pub fn createConvolutionMatrix(
    allocator: Allocator,
    kernel: []const f64,
    kernel_radius: usize,
    signal_size: usize,
) !SymmetricBandMatrix {
    var matrix = try SymmetricBandMatrix.init(allocator, signal_size, kernel_radius);
    matrix.clear();

    // Construct band matrix from symmetric kernel
    for (0..signal_size) |i| {
        for (0..@min(kernel_radius + 1, signal_size - i)) |offset| {
            const j = i + offset;
            if (j < signal_size and offset < kernel.len) {
                matrix.set(i, j, kernel[offset]);
            }
        }
    }

    return matrix;
}

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
    try eqa(2.0, matrix.get(0, 0), 1e-10);
    try eqa(1.0, matrix.get(0, 1), 1e-10);
    try eqa(1.0, matrix.get(1, 0), 1e-10); // symmetric
    try eqa(2.0, matrix.get(1, 1), 1e-10);
    try eqa(1.0, matrix.get(1, 2), 1e-10);
    try eqa(1.0, matrix.get(2, 1), 1e-10); // symmetric
    try eqa(2.0, matrix.get(2, 2), 1e-10);

    // Test out-of-band elements
    try eqa(0.0, matrix.get(0, 2), 1e-10);
    try eqa(0.0, matrix.get(2, 0), 1e-10);
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

test "convolution matrix creation" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    const kernel = [_]f64{ 0.25, 0.5, 0.25 };
    const signal_size = 5;

    var matrix = try createConvolutionMatrix(allocator, &kernel, 2, signal_size);
    defer matrix.deinit(allocator);

    // Check that the kernel values are properly stored
    try eqa(0.25, matrix.get(0, 0), 1e-10);
    try eqa(0.5, matrix.get(0, 1), 1e-10);
    try eqa(0.25, matrix.get(0, 2), 1e-10);

    // Test symmetry
    try eqa(0.5, matrix.get(1, 0), 1e-10);
    try eqa(0.25, matrix.get(1, 1), 1e-10);
    try eqa(0.5, matrix.get(1, 2), 1e-10);
}

test "matrix clear operation" {
    const allocator = testing.allocator;
    const eqa = testing.expectApproxEqAbs;

    var matrix = try SymmetricBandMatrix.init(allocator, 3, 1);
    defer matrix.deinit(allocator);

    // Set some values
    matrix.set(0, 0, 5.0);
    matrix.set(0, 1, 3.0);

    // Clear the matrix
    matrix.clear();

    // Check that all values are zero
    try eqa(0.0, matrix.get(0, 0), 1e-10);
    try eqa(0.0, matrix.get(0, 1), 1e-10);
    try eqa(0.0, matrix.get(1, 0), 1e-10);
}
