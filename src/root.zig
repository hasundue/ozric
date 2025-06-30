//! By convention, root.zig is the root source file when making a library. If
//! you are making an executable, the convention is to delete this file and
//! start with main.zig instead.
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
