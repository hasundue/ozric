const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const ceres = b.dependency("ceres", .{
        .target = target,
        .optimize = optimize,
    });

    const eigen = b.dependency("eigen", .{
        .target = target,
        .optimize = optimize,
    });

    const lib_mod = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "ozric",
        .root_module = lib_mod,
    });

    lib.addCSourceFiles(.{
        .root = b.path("src"),
        .files = &.{"solver.cc"},
    });

    // Add all ceres source files recursively with MINIGLOG for 2.2.0
    const ceres_flags = &.{
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

    // Add ceres source files with required threading support
    const ceres_sources = [_][]const u8{
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

    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres"),
        .files = &ceres_sources,
        .flags = ceres_flags,
    });

    // Since we set CERES_RESTRICT_SCHUR_SPECIALIZATION, we only need to
    // add the dynamic-sized Schur eliminator and partitioned matrix view
    // template specializations.
    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres/generated"),
        .files = &.{
            "schur_eliminator_d_d_d.cc",
            "partitioned_matrix_view_d_d_d.cc",
        },
        .flags = ceres_flags,
    });

    // Add miniglog source file with proper export macros
    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres/miniglog/glog"),
        .files = &.{"logging.cc"},
        .flags = &.{
            "-std=c++17",
            "-DCERES_EXPORT=", // Define CERES_EXPORT as empty
            "-DCERES_NO_EXPORT=", // Define CERES_NO_EXPORT as empty
            "-DMINIGLOG", // Use minimal logging
        },
    });

    // Add ceres include directories
    lib.addIncludePath(ceres.path("include"));
    lib.addIncludePath(ceres.path(".")); // For internal headers
    lib.addIncludePath(ceres.path("internal")); // For internal ceres headers
    lib.addIncludePath(ceres.path("config")); // For config headers like export.h
    lib.addIncludePath(ceres.path("internal/ceres/miniglog")); // For miniglog headers

    // Add Eigen include directory
    lib.addIncludePath(eigen.path(".")); // Eigen is header-only

    // Add our local include directory
    lib.addIncludePath(b.path("include"));

    // MINIGLOG eliminates the need for Abseil dependencies
    // lib.addIncludePath(abseil.path("."));

    lib.linkLibCpp();

    b.installArtifact(lib);

    const exe_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    exe_mod.addImport("ozric_lib", lib_mod);

    const exe = b.addExecutable(.{
        .name = "ozric",
        .root_module = exe_mod,
    });

    // We don't need to add C++ source files to the executable since
    // the library already includes them.

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    const lib_unit_tests = b.addTest(.{
        .root_module = lib_mod,
    });

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const exe_unit_tests = b.addTest(.{
        .root_module = exe_mod,
    });

    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
    test_step.dependOn(&run_exe_unit_tests.step);

    const wasm_step = b.step("wasm", "Build WASM using Emscripten");

    // Create Emscripten compilation command
    const emcc_cmd = b.addSystemCommand(&.{
        "emcc",
        "-std=c++17",
        "-pthread",
        "-sASYNCIFY",
        "-sPTHREAD_POOL_SIZE=4",
        "-sINITIAL_MEMORY=167772160", // 160MB
        "-sUSE_OFFSET_CONVERTER",
        "-sEXPORTED_FUNCTIONS=_test_ceres,_run_hello_world",
        "-sEXPORTED_RUNTIME_METHODS=ccall,cwrap,PThread",
        "-O0",
        "--cache",
        ".zig-cache/emcc",
        // Ceres defines
        "-DCERES_NO_SUITESPARSE",
        "-DCERES_NO_CXSPARSE",
        "-DCERES_EXPORT=",
        "-DCERES_NO_EXPORT=",
        "-DCERES_NO_PROTOCOL_BUFFERS",
        "-DCERES_NO_CUDA",
        "-DCERES_NO_LAPACK",
        "-DCERES_METIS_VERSION=\"5.1.0\"",
        "-DMINIGLOG",
        "-DCERES_RESTRICT_SCHUR_SPECIALIZATION",
        "-DCERES_NO_ACCELERATE_SPARSE",
        // Include paths
    });

    // Add include paths
    emcc_cmd.addArg(b.fmt("-I{s}", .{ceres.path("include").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{ceres.path(".").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{ceres.path("internal").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{ceres.path("config").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{ceres.path("internal/ceres/miniglog").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{eigen.path(".").getPath(b)}));
    emcc_cmd.addArg(b.fmt("-I{s}", .{b.path("include").getPath(b)}));

    // Add our C++ wrapper
    emcc_cmd.addArg(b.fmt("{s}/src/solver.cc", .{b.build_root.path.?}));

    // Add all essential Ceres sources
    for (ceres_sources) |source| {
        emcc_cmd.addArg(b.fmt("{s}/internal/ceres/{s}", .{ ceres.path(".").getPath(b), source }));
    }

    // Add template specializations
    emcc_cmd.addArg(b.fmt("{s}/internal/ceres/generated/schur_eliminator_d_d_d.cc", .{ceres.path(".").getPath(b)}));
    emcc_cmd.addArg(b.fmt("{s}/internal/ceres/generated/partitioned_matrix_view_d_d_d.cc", .{ceres.path(".").getPath(b)}));

    // Add miniglog
    emcc_cmd.addArg(b.fmt("{s}/internal/ceres/miniglog/glog/logging.cc", .{ceres.path(".").getPath(b)}));

    // Output file
    emcc_cmd.addArg("-o");
    emcc_cmd.addArg("zig-out/bin/ozric_wasm.js");

    wasm_step.dependOn(&emcc_cmd.step);
}
