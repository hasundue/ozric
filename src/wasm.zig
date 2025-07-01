const std = @import("std");
const lib = @import("ozric_lib");

export fn getCeresVersion() i32 {
    return lib.getCeresVersion();
}

export fn runHelloWorld() f64 {
    return lib.runHelloWorld();
}

// Main function for WASM executable - test what we can
pub fn main() !void {
    std.log.info("Ozric WASM module loaded", .{});

    // Test basic functionality
    const version = getCeresVersion();
    std.log.info("Ceres version: {}", .{version});

    // For now, just test if the function exists (may fail due to missing deps)
    std.log.info("Attempting to run hello world...", .{});
    // Commented out until we have all dependencies
    // const result = runHelloWorld();
    // std.log.info("Hello world result: {d} (should be ~10.0)", .{result});
}
