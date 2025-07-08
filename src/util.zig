const std = @import("std");
const t = std.testing;

/// Returns the next power of two greater than or equal to n
/// Used for determining FFT buffer sizes
pub fn nextPowerOfTwo(n: usize) usize {
    if (n == 0) return 1;

    var power: usize = 1;
    while (power < n) {
        power <<= 1;
    }
    return power;
}

test nextPowerOfTwo {
    try t.expectEqual(@as(usize, 1), nextPowerOfTwo(0));
    try t.expectEqual(@as(usize, 1), nextPowerOfTwo(1));
    try t.expectEqual(@as(usize, 4), nextPowerOfTwo(3));
    try t.expectEqual(@as(usize, 8), nextPowerOfTwo(5));
    try t.expectEqual(@as(usize, 16), nextPowerOfTwo(16));
}
