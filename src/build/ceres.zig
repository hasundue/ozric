const std = @import("std");

/// Compilation flags for Ceres solver with MINIGLOG for 2.2.0
pub const ceres_flags = [_][]const u8{
    "-std=c++17",
    "-DCERES_NO_SUITESPARSE",
    "-DCERES_NO_CXSPARSE",
    "-DCERES_EXPORT=", // Define CERES_EXPORT as empty
    "-DCERES_NO_EXPORT=", // Define CERES_NO_EXPORT as empty
    "-DCERES_NO_PROTOCOL_BUFFERS",
    "-DCERES_NO_THREADS",
    "-DCERES_NO_CUDA", // Disable CUDA support
    "-DCERES_NO_LAPACK", // Disable LAPACK dependencies
    "-DCERES_METIS_VERSION=\"5.1.0\"", // Define METIS version
    "-DMINIGLOG", // Use minimal logging instead of full glog
    "-DCERES_RESTRICT_SCHUR_SPECIALIZATION", // Restrict Schur specializations
    "-DCERES_NO_ACCELERATE_SPARSE", // Disable Accelerate sparse support
};

/// Core Ceres source files in internal/ceres/
const ceres_core_sources = [_][]const u8{
    // Core utilities
    "array_utils.cc",
    "stringprintf.cc",
    "wall_time.cc",

    // Problem core
    "problem.cc",
    "problem_impl.cc",
    "cost_function.cc",
    "loss_function.cc",
    "corrector.cc",
    "residual_block.cc",

    // Solver core and utilities
    "solver.cc",
    "minimizer.cc",
    "trust_region_minimizer.cc",
    "trust_region_strategy.cc",
    "levenberg_marquardt_strategy.cc",
    "types.cc",

    // Basic linear solver (QR only)
    "linear_solver.cc",
    "dense_qr_solver.cc",
    "dense_qr.cc",

    // Program structure and utilities
    "program.cc",
    "preprocessor.cc",
    "manifold.cc",
    "detect_structure.cc",

    // Evaluation basics
    "evaluator.cc",
    "block_evaluate_preparer.cc",
    "block_jacobian_writer.cc",
    "compressed_row_jacobian_writer.cc",
    "block_sparse_matrix.cc",
    "compressed_row_sparse_matrix.cc",
    "triplet_sparse_matrix.cc",

    // Threading support (required)
    "thread_pool.cc",

    // Sparse matrix base class and linear operators
    "sparse_matrix.cc",
    "linear_operator.cc",

    // Mathematical utilities
    "gradient_checker.cc",
    "is_close.cc",

    // Dense matrix support
    "dense_sparse_matrix.cc",

    // Parallel operations
    "parallel_vector_ops.cc",
    "parallel_invoke.cc",
    "parallel_utils.cc",

    // Solver utilities and templates
    "solver_utils.cc",
    "schur_templates.cc",

    // Block structure and jacobian support
    "block_structure.cc",
    "scratch_evaluate_preparer.cc",
    "dynamic_compressed_row_jacobian_writer.cc",

    // Parameter ordering and preprocessing
    "parameter_block_ordering.cc",
    "trust_region_preprocessor.cc",
    "reorder_program.cc",
    "coordinate_descent_minimizer.cc",

    // Callback support
    "iteration_callback.cc",
    "gradient_problem_solver.cc",
    "callbacks.cc",
    "evaluation_callback.cc",

    // Gradient problem support
    "gradient_problem.cc",

    // Linear least squares support
    "linear_least_squares_problems.cc",

    // Dense matrix linear algebra
    "dense_cholesky.cc",
    "dense_normal_cholesky_solver.cc",

    // Sparse matrix operations
    "dynamic_compressed_row_sparse_matrix.cc",
    "sparse_cholesky.cc",

    // Preconditioner systems
    "preconditioner.cc",
    "block_jacobi_preconditioner.cc",
    "subset_preconditioner.cc",
    "schur_jacobi_preconditioner.cc",
    "power_series_expansion_preconditioner.cc",

    // Additional linear solvers
    "cgnr_solver.cc",
    "iterative_schur_complement_solver.cc",
    "sparse_normal_cholesky_solver.cc",
    "dynamic_sparse_normal_cholesky_solver.cc",

    // Line search and optimization
    "line_search.cc",
    "line_search_direction.cc",
    "line_search_minimizer.cc",
    "line_search_preprocessor.cc",

    // Advanced mathematics
    "polynomial.cc",
    "function_sample.cc",
    "trust_region_step_evaluator.cc",
    "low_rank_inverse_hessian.cc",

    // Schur complement methods
    "schur_eliminator.cc",
    "implicit_schur_complement.cc",

    // Block matrix operations
    "block_random_access_diagonal_matrix.cc",
    "block_random_access_sparse_matrix.cc",
    "block_random_access_dense_matrix.cc",

    // Computational utilities
    "inner_product_computer.cc",
    "iterative_refiner.cc",
    "compressed_col_sparse_matrix_utils.cc",
    "residual_block_utils.cc",

    // File I/O operations
    "file.cc",

    // Gradient checking support
    "gradient_checking_cost_function.cc",

    // Context
    "context.cc",
    "context_impl.cc",

    // Additional algorithm implementations
    "dogleg_strategy.cc",
    "schur_complement_solver.cc",
    "visibility_based_preconditioner.cc",

    // Visibility and clustering algorithms
    "visibility.cc",
    "canonical_views_clustering.cc",
    "single_linkage_clustering.cc",
    "block_random_access_matrix.cc",

    // Template specialization base classes
    "partitioned_matrix_view.cc",
};

/// Template specialization sources in internal/ceres/generated/
const ceres_generated_sources = [_][]const u8{
    // Template specializations (dynamic-sized only due to CERES_RESTRICT_SCHUR_SPECIALIZATION)
    "schur_eliminator_d_d_d.cc",
    "partitioned_matrix_view_d_d_d.cc",
};

/// Miniglog sources in internal/ceres/miniglog/glog/
const miniglog_sources = [_][]const u8{
    "logging.cc",
};

/// Comprehensive list of all Ceres source files with full paths
pub const all_ceres_sources = blk: {
    var sources: [ceres_core_sources.len + ceres_generated_sources.len + miniglog_sources.len][]const u8 = undefined;
    var i: usize = 0;

    // Add core sources with internal/ceres/ prefix
    for (ceres_core_sources) |source| {
        sources[i] = "internal/ceres/" ++ source;
        i += 1;
    }

    // Add generated sources with internal/ceres/generated/ prefix
    for (ceres_generated_sources) |source| {
        sources[i] = "internal/ceres/generated/" ++ source;
        i += 1;
    }

    // Add miniglog sources with internal/ceres/miniglog/glog/ prefix
    for (miniglog_sources) |source| {
        sources[i] = "internal/ceres/miniglog/glog/" ++ source;
        i += 1;
    }

    break :blk sources;
};

/// Add Ceres-solver support to a library or executable
pub fn addCeresSupport(b: *std.Build, artifact: *std.Build.Step.Compile, ceres: *std.Build.Dependency, eigen: *std.Build.Dependency) void {
    // Add all ceres source files (core, generated, and miniglog) with unified flags
    artifact.addCSourceFiles(.{
        .root = ceres.path("."),
        .files = &all_ceres_sources,
        .flags = &ceres_flags,
    });

    // Add ceres include directories
    artifact.addIncludePath(ceres.path("include"));
    artifact.addIncludePath(ceres.path(".")); // For internal headers
    artifact.addIncludePath(ceres.path("internal")); // For internal ceres headers
    artifact.addIncludePath(ceres.path("config")); // For config headers like export.h
    artifact.addIncludePath(ceres.path("internal/ceres/miniglog")); // For miniglog headers

    // Add Eigen include directory
    artifact.addIncludePath(eigen.path(".")); // Eigen is header-only

    // Add our local include directory
    artifact.addIncludePath(b.path("include"));

    // Link C++ standard library
    artifact.linkLibCpp();
}

pub const WasmSourceInfo = struct {
    path: []const u8,
    base_dir: std.Build.LazyPath,
};

/// Get all Ceres source files for WASM compilation
pub fn getCeresSources(ceres: *std.Build.Dependency, allocator: std.mem.Allocator) !std.ArrayList(WasmSourceInfo) {
    var sources = std.ArrayList(WasmSourceInfo).init(allocator);

    // Add all ceres sources (core, generated, and miniglog)
    for (all_ceres_sources) |source| {
        try sources.append(.{
            .path = source,
            .base_dir = ceres.path("."),
        });
    }

    return sources;
}

/// Compilation command for Emscripten C++ files with Ceres support
pub const emcc_compile_command = [_][]const u8{
    "ccache",
    "emcc",
    "-c",
    "-pthread",
    "-O0",
    "--cache",
    ".zig-cache/emcc",
} ++ ceres_flags;

/// Add emcc include paths to a system command
pub fn addEmccIncludePaths(compile_cmd: *std.Build.Step.Run, b: *std.Build, ceres: *std.Build.Dependency, eigen: *std.Build.Dependency) void {
    // Allowlist absolute paths to suppress warnings
    compile_cmd.addArg(b.fmt("--valid-abspath={s}", .{ceres.path(".").getPath(b)}));
    compile_cmd.addArg(b.fmt("--valid-abspath={s}", .{eigen.path(".").getPath(b)}));
    compile_cmd.addArg(b.fmt("--valid-abspath={s}", .{b.path("include").getPath(b)}));

    // Add include paths
    compile_cmd.addArg(b.fmt("-I{s}", .{ceres.path("include").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{ceres.path(".").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{ceres.path("internal").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{ceres.path("config").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{ceres.path("internal/ceres/miniglog").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{eigen.path(".").getPath(b)}));
    compile_cmd.addArg(b.fmt("-I{s}", .{b.path("include").getPath(b)}));
}

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
