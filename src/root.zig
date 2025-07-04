const std = @import("std");
const testing = std.testing;

const c = @cImport({
    @cInclude("solver.h");
});

pub fn getCeresVersion() i32 {
    return c.test_ceres();
}

pub fn runHelloWorld() f64 {
    return c.run_hello_world();
}

pub export fn add(a: i32, b: i32) i32 {
    return a + b;
}

test "basic add functionality" {
    try testing.expect(add(3, 7) == 10);
}

test "ozric main functionality" {
    const ceres_version = getCeresVersion();

    std.debug.print("Ozric - Ornstein-Zernike equation solver\n", .{});
    std.debug.print("Using Ceres Solver version: {d}\n", .{ceres_version});

    const result = runHelloWorld();
    std.debug.print("Hello World result: x = {d} (should be ~10.0)\n", .{result});

    try testing.expectApproxEqAbs(10.0, result, 1e-4);
}
