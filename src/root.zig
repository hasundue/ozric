pub const convolution = @import("convolution.zig");

test "lib" {
    _ = @import("grid.zig");
    _ = @import("hardsphere.zig");
    _ = @import("convolution.zig");
}

test "linear" {
    _ = @import("linear/dsbmv.zig");
    _ = @import("linear/sb.zig");
}
