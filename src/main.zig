const std = @import("std");
const ozric = @import("ozric");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.log.info("Ozric - Ornstein-Zernike Equation Solver", .{});
    
    // Example: Initialize a simple solver
    var solver = ozric.Solver.init(allocator);
    defer solver.deinit();
    
    std.log.info("Solver initialized successfully", .{});
}