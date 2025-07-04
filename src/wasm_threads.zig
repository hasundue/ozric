const std = @import("std");
const lib = @import("ozric_lib");

export fn getCeresVersion() i32 {
    return lib.getCeresVersion();
}

export fn runHelloWorld() i32 {
    return lib.runHelloWorld();
}

export fn testOptimization() i32 {
    return lib.testOptimization();
}

pub fn main() !void {
    std.log.info("Ozric WASM Threads module loaded with full Ceres support", .{});

    const version = getCeresVersion();
    std.log.info("Ceres version: {}", .{version});

    std.log.info("Running full hello world optimization...", .{});
    const hello_result = runHelloWorld();
    std.log.info("Hello world result: {}", .{hello_result});

    std.log.info("Running optimization test...", .{});
    const opt_result = testOptimization();
    std.log.info("Optimization test result: {}", .{opt_result});
}
