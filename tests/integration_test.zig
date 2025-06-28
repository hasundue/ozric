const std = @import("std");
const ozric = @import("ozric");

test "integration: basic solver workflow" {
    const allocator = std.testing.allocator;

    var solver = ozric.Solver.init(allocator);
    defer solver.deinit();

    // Test with typical liquid parameters
    const density = 0.8; // g/cm³
    const temperature = 298.15; // K

    try solver.solve(density, temperature);
}

test "integration: vector operations with BLAS" {
    const allocator = std.testing.allocator;

    var v1 = try ozric.Vector.init(allocator, 100);
    defer v1.deinit();

    var v2 = try ozric.Vector.init(allocator, 100);
    defer v2.deinit();

    // Initialize with some test data
    for (0..100) |i| {
        v1.data[i] = @floatFromInt(i + 1);
        v2.data[i] = 1.0;
    }

    const result = v1.dot(&v2);
    try std.testing.expectEqual(@as(f64, 5050.0), result); // sum of 1 to 100

    // Test scaling
    v1.scale(2.0);
    try std.testing.expectEqual(@as(f64, 2.0), v1.data[0]);
    try std.testing.expectEqual(@as(f64, 200.0), v1.data[99]);
}
