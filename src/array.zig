const std = @import("std");

pub fn allocArraySeries(
    comptime T: type,
    allocator: std.mem.Allocator,
    n: usize,
    size: usize,
) ![][]T {
    const data = try allocator.alloc(T, n * size);

    const row_ptrs = try allocator.alloc([]T, n);

    for (row_ptrs, 0..) |*row_ptr, i| {
        const start = i * size;
        row_ptr.* = data[start .. start + size];
    }

    return row_ptrs;
}

pub fn freeArraySeries(comptime T: type, allocator: std.mem.Allocator, array: [][]T) void {
    if (array.len > 0 and array[0].len > 0) {
        const total_size = array.len * array[0].len;
        const data_ptr = array[0].ptr;
        allocator.free(data_ptr[0..total_size]);
    }
    allocator.free(array);
}

test "array series allocation and deallocation" {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const array = try allocArraySeries(f64, allocator, 3, 4);
    defer freeArraySeries(f64, allocator, array);

    try std.testing.expect(array.len == 3);
    try std.testing.expect(array[0].len == 4);
    try std.testing.expect(array[1].len == 4);
    try std.testing.expect(array[2].len == 4);

    array[0][0] = 1.0;
    array[2][3] = 2.0;
    try std.testing.expect(array[0][0] == 1.0);
    try std.testing.expect(array[2][3] == 2.0);
}
