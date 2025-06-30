const std = @import("std");

/// Programmatically generates the list of Ceres template specialization source files.
///
/// This function creates filenames for the generated template specializations that
/// Ceres uses for optimized Schur complement elimination. The specializations cover
/// various combinations of block sizes commonly found in bundle adjustment problems.
///
/// Template parameters: (row_block_size, e_block_size, f_block_size)
/// - 2,3,4: Fixed small block sizes for efficient compile-time optimization
/// - d: Dynamic block size for flexibility with variable-sized problems
/// - 8,9: Additional sizes for specific computer vision applications
///
/// Returns a slice of filenames that must be freed by the caller.
fn generateTemplateSpecializations(allocator: std.mem.Allocator) ![][]const u8 {
    var sources = std.ArrayList([]const u8).init(allocator);

    // Template parameter combinations used by Ceres for bundle adjustment optimization
    // These correspond to common block structures in computer vision problems
    const template_params = [_][]const u8{
        "2_2_2", "2_2_3", "2_2_4", "2_2_d", // 2x2 camera blocks
        "2_3_3", "2_3_4", "2_3_9", "2_3_d", // 2x3 camera blocks
        "2_4_3", "2_4_4", "2_4_8", "2_4_9", "2_4_d", // 2x4 camera blocks
        "2_d_d", "4_4_2", "4_4_3", "4_4_4", "4_4_d", // Dynamic and 4x4 blocks
        "d_d_d", // Fully dynamic
    };

    // Generate SchurEliminator template specializations
    for (template_params) |params| {
        const filename = try std.fmt.allocPrint(allocator, "schur_eliminator_{s}.cc", .{params});
        try sources.append(filename);
    }

    // Generate PartitionedMatrixView template specializations
    for (template_params) |params| {
        const filename = try std.fmt.allocPrint(allocator, "partitioned_matrix_view_{s}.cc", .{params});
        try sources.append(filename);
    }

    return sources.toOwnedSlice();
}
