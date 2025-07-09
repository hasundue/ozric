test "lib" {
    _ = @import("grid.zig");
    _ = @import("hardsphere.zig");
    _ = @import("convolution.zig");
}

test "la" {
    _ = @import("la/dsbmv.zig");
}
