const std = @import("std");
const lib = @import("ozric_lib");

export fn getCeresVersion() i32 {
    return lib.getCeresVersion();
}

export fn runHelloWorld() f64 {
    return lib.runHelloWorld();
}

// Main function for WASM executable
pub fn main() !void {
    std.log.info("Ozric WASM module loaded", .{});
}
