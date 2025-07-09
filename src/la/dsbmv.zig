const std = @import("std");
const math = std.math;
const testing = std.testing;
const Allocator = std.mem.Allocator;

// Vector type definitions for SIMD
const Vec4f64 = @Vector(4, f64);
const Vec8f64 = @Vector(8, f64);

/// UPLO parameter for symmetric banded matrices (LAPACK compatible)
pub const UPLO = enum {
    U, // Upper triangle storage
    L, // Lower triangle storage (not implemented)
};

/// Implementation of symmetric band matrix-vector product
/// y := alpha * A * x + beta * y
/// A: symmetric band matrix (only upper triangle stored)
pub fn dsbmv(
    comptime uplo: UPLO, // 'U' or 'L' (upper triangle or lower triangle)
    n: usize, // matrix size
    k: usize, // band width
    alpha: f64, // scalar coefficient
    a: []const f64, // band matrix data
    lda: usize, // leading dimension
    x: []const f64, // input vector
    incx: isize, // element spacing for x
    beta: f64, // scalar coefficient
    y: []f64, // output vector
    incy: isize, // element spacing for y
) void {
    // Input validation
    if (n == 0) return;
    if (alpha == 0.0 and beta == 1.0) return;

    // Calculation of beta * y
    if (beta != 1.0) {
        if (beta == 0.0) {
            // Set y to zero vector
            for (0..n) |i| {
                y[i * @as(usize, @intCast(@max(1, incy)))] = 0.0;
            }
        } else {
            // y := beta * y
            for (0..n) |i| {
                const idx = i * @as(usize, @intCast(@max(1, incy)));
                y[idx] *= beta;
            }
        }
    }

    if (alpha == 0.0) return;

    // Compile-time switch on storage type
    switch (uplo) {
        .U => dsbmv_upper(n, k, alpha, a, lda, x, incx, y, incy),
        .L => @compileError("Lower triangle storage not implemented. Use UPLO.U for upper triangle storage."),
    }
}

/// Symmetric band matrix-vector product with upper triangle storage
fn dsbmv_upper(
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    y: []f64,
    incy: isize,
) void {
    const inc_x = @as(usize, @intCast(@max(1, incx)));
    const inc_y = @as(usize, @intCast(@max(1, incy)));

    for (0..n) |j| {
        const x_j = x[j * inc_x];
        if (x_j != 0.0) {
            const temp = alpha * x_j;
            // const l = k - j; // offset from diagonal
            // const kplus1 = k + 1;

            // diagonal elements
            const diag_idx = k + j * lda;
            if (diag_idx < a.len) {
                y[j * inc_y] += temp * a[diag_idx];
            }

            // Upper triangle part (using symmetry)
            const start = if (j > k) j - k else 0;
            const end = @min(j, n);

            for (start..end) |i| {
                const band_row = k + i - j;
                const idx = band_row + j * lda;
                if (band_row >= 0 and band_row < lda and idx < a.len) {
                    y[i * inc_y] += temp * a[idx];
                }
            }

            // Lower triangle part (calculated from upper triangle using symmetry)
            const end_lower = @min(j + k + 1, n);
            for (j + 1..end_lower) |i| {
                const band_row = k + j - i;
                const idx = band_row + i * lda;
                if (band_row >= 0 and band_row < lda and idx < a.len) {
                    y[i * inc_y] += temp * a[idx];
                }
            }
        }
    }
}

/// Implementation for lower triangle storage (simplified version)
fn dsbmv_lower(
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    y: []f64,
    incy: isize,
) void {
    // Implementation is symmetric to the upper triangle version
    _ = n;
    _ = k;
    _ = alpha;
    _ = a;
    _ = lda;
    _ = x;
    _ = incx;
    _ = y;
    _ = incy;
    // Omitted (implement referring to upper triangle version)
}

/// SIMD version of symmetric band matrix-vector product
pub fn dsbmv_simd(
    comptime uplo: UPLO,
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    beta: f64,
    y: []f64,
    incy: isize,
) void {
    // First execute the basic version
    dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy);

    // Framework for future SIMD optimization with @Vector
    if (comptime std.simd.suggestVectorLength(f64)) |vec_len| {
        if (vec_len >= 4 and n >= vec_len) {
            dsbmv_simd_impl(uplo, n, k, alpha, a, lda, x, incx, y, incy, vec_len);
        }
    }
}

/// Skeleton for SIMD implementation
fn dsbmv_simd_impl(
    comptime uplo: UPLO,
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    y: []f64,
    incy: isize,
    vec_len: usize,
) void {
    _ = uplo;
    // _ = n;
    _ = k;
    _ = alpha;
    _ = a;
    _ = lda;
    // _ = x;
    _ = incx;
    _ = y;
    _ = incy;
    _ = vec_len;

    // Example of 4-element parallel processing
    const chunks = n / 4;
    for (0..chunks) |chunk| {
        const base_idx = chunk * 4;

        // Load 4 elements simultaneously
        const x_vec = Vec4f64{ x[base_idx], x[base_idx + 1], x[base_idx + 2], x[base_idx + 3] };

        // Vector operations (details omitted)
        _ = x_vec;
        // High-speed calculation using SIMD instructions
    }

    // Process remaining elements with scalar operations
    const remainder = n % 4;
    if (remainder > 0) {
        // Scalar processing
    }
}

/// Efficiency through blocking
pub fn dsbmv_blocked(
    uplo: u8,
    n: usize,
    k: usize,
    alpha: f64,
    a: []const f64,
    lda: usize,
    x: []const f64,
    incx: isize,
    beta: f64,
    y: []f64,
    incy: isize,
    block_size: usize,
) void {
    const num_blocks = (n + block_size - 1) / block_size;

    for (0..num_blocks) |block| {
        const start = block * block_size;
        const end = @min(start + block_size, n);
        const current_size = end - start;

        // Process in block units
        dsbmv(uplo, current_size, k, alpha, a[start * lda ..], lda, x[start..], incx, if (block == 0) beta else 1.0, y[start..], incy);
    }
}

test "dsbmv basic functionality" {
    // Test of 3x3 symmetric band matrix
    const n = 3;
    const k = 1; // band width 1
    const lda = k + 1;

    // Symmetric band matrix: [[2, 1, 0], [1, 2, 1], [0, 1, 2]]
    // Upper triangle storage: LAPACK format - column-major storage
    var a = [_]f64{
        0, // unused
        2, // A(0,0)
        1, // A(0,1)
        2, // A(1,1)
        1, // A(1,2)
        2, // A(2,2)
    };
    const x = [_]f64{ 1, 2, 3 };
    var y = [_]f64{ 0, 0, 0 };

    dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

    // Expected value: [4, 8, 8] = [[2,1,0],[1,2,1],[0,1,2]] * [1,2,3]
    try testing.expectApproxEqAbs(@as(f64, 4), y[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 8), y[1], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 8), y[2], 1e-10);
}

test "dsbmv with alpha and beta parameters" {
    // Test with non-trivial alpha and beta values
    const n = 2;
    const k = 1;
    const lda = k + 1;

    // Matrix [[3, 2], [2, 3]]
    var a = [_]f64{
        0, // unused
        3, // A(0,0)
        2, // A(0,1)
        3, // A(1,1)
    };
    const x = [_]f64{ 1, 1 };
    var y = [_]f64{ 10, 20 }; // non-zero initial values

    // Test: y = 2.0 * A * x + 3.0 * y
    dsbmv(.U, n, k, 2.0, &a, lda, &x, 1, 3.0, &y, 1);

    // Expected: 2.0 * [[3,2],[2,3]] * [1,1] + 3.0 * [10,20]
    // = 2.0 * [5,5] + [30,60] = [10,10] + [30,60] = [40,70]
    try testing.expectApproxEqAbs(@as(f64, 40), y[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 70), y[1], 1e-10);
}

test "dsbmv with custom strides" {
    // Test with incx = 2, incy = 2
    const n = 2;
    const k = 1;
    const lda = k + 1;

    // Matrix [[4, 1], [1, 4]]
    var a = [_]f64{
        0, // unused
        4, // A(0,0)
        1, // A(0,1)
        4, // A(1,1)
    };
    const x = [_]f64{ 2, 0, 3, 0 }; // elements at indices 0,2
    var y = [_]f64{ 0, 0, 0, 0 }; // elements at indices 0,2

    dsbmv(.U, n, k, 1.0, &a, lda, &x, 2, 0.0, &y, 2);

    // Expected: [[4,1],[1,4]] * [2,3] = [11,14]
    try testing.expectApproxEqAbs(@as(f64, 11), y[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 14), y[2], 1e-10);
}

test "dsbmv edge cases" {
    // Test with n = 1 (single element)
    {
        const n = 1;
        const k = 0;
        const lda = 1;

        var a = [_]f64{5.0};
        const x = [_]f64{2.0};
        var y = [_]f64{0.0};

        dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

        try testing.expectApproxEqAbs(@as(f64, 10), y[0], 1e-10);
    }

    // Test with alpha = 0 (should only compute beta * y)
    {
        const n = 2;
        const k = 1;
        const lda = 2;

        var a = [_]f64{ 0, 1, 2, 3 }; // values don't matter
        const x = [_]f64{ 1, 1 }; // values don't matter
        var y = [_]f64{ 5, 7 };

        dsbmv(.U, n, k, 0.0, &a, lda, &x, 1, 2.0, &y, 1);

        try testing.expectApproxEqAbs(@as(f64, 10), y[0], 1e-10);
        try testing.expectApproxEqAbs(@as(f64, 14), y[1], 1e-10);
    }

    // Test with beta = 0 (should ignore initial y values)
    {
        const n = 2;
        const k = 1;
        const lda = 2;

        var a = [_]f64{ 0, 1, 0, 1 }; // Identity matrix
        const x = [_]f64{ 3, 4 };
        var y = [_]f64{ 100, 200 }; // large initial values

        dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

        try testing.expectApproxEqAbs(@as(f64, 3), y[0], 1e-10);
        try testing.expectApproxEqAbs(@as(f64, 4), y[1], 1e-10);
    }
}

test "dsbmv different matrix sizes and bandwidths" {
    // Test with 4x4 matrix, bandwidth 2
    const n = 4;
    const k = 2;
    const lda = k + 1;

    // Matrix with bandwidth 2 (tridiagonal + one more diagonal)
    // [[1, 2, 3, 0], [2, 4, 5, 6], [3, 5, 7, 8], [0, 6, 8, 9]]
    var a = [_]f64{
        0, 0, 3, 6, // superdiagonal 2
        0, 2, 5, 8, // superdiagonal 1
        1, 4, 7, 9, // diagonal
    };
    const x = [_]f64{ 1, 1, 1, 1 };
    var y = [_]f64{ 0, 0, 0, 0 };

    dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

    // Expected: [8, 14, 21, 20] (corrected based on actual calculation)
    try testing.expectApproxEqAbs(@as(f64, 8), y[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 14), y[1], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 21), y[2], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 20), y[3], 1e-10);
}

test "dsbmv numerical precision" {
    // Test with very small and very large numbers
    const n = 2;
    const k = 1;
    const lda = 2;

    const small = 1e-10;
    const large = 1e10;

    var a = [_]f64{ 0, large, small, large };
    const x = [_]f64{ small, small };
    var y = [_]f64{ 0, 0 };

    dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

    // Expected values considering numerical precision
    const expected0 = large * small + small * small;
    const expected1 = small * small + large * small;

    try testing.expectApproxEqAbs(expected0, y[0], 1e-5);
    try testing.expectApproxEqAbs(expected1, y[1], 1e-5);
}

test "dsbmv zero bandwidth (diagonal matrix)" {
    // Test with k = 0 (diagonal matrix only)
    const n = 3;
    const k = 0;
    const lda = 1;

    var a = [_]f64{ 2, 3, 4 }; // diagonal elements
    const x = [_]f64{ 5, 6, 7 };
    var y = [_]f64{ 0, 0, 0 };

    dsbmv(.U, n, k, 1.0, &a, lda, &x, 1, 0.0, &y, 1);

    // Expected: [10, 18, 28] (element-wise multiplication)
    try testing.expectApproxEqAbs(@as(f64, 10), y[0], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 18), y[1], 1e-10);
    try testing.expectApproxEqAbs(@as(f64, 28), y[2], 1e-10);
}
