const std = @import("std");
const testing = std.testing;

const c = @cImport({
    @cInclude("solver.h");
});

test "solver" {
    const ceres_version = c.test_ceres();

    std.debug.print("Ozric - Ornstein-Zernike equation solver\n", .{});
    std.debug.print("Using Ceres Solver version: {d}\n", .{ceres_version});

    const result = c.solve();
    std.debug.print("Solver result: x = {d} (should be ~10.0)\n", .{result});

    try testing.expectApproxEqAbs(10.0, result, 1e-4);
}
