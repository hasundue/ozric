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
    const essential_ceres_sources = [_][]const u8{
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
        .files = &essential_ceres_sources,
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

    // Create WASM executable with minimal library module (no library artifact)
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .wasi,
    });

    const wasm_lib_mod = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = wasm_target,
        .optimize = optimize,
    });

    // Add include paths to the WASM library module for C imports
    wasm_lib_mod.addIncludePath(ceres.path("include"));
    wasm_lib_mod.addIncludePath(ceres.path("."));
    wasm_lib_mod.addIncludePath(ceres.path("internal"));
    wasm_lib_mod.addIncludePath(ceres.path("config"));
    wasm_lib_mod.addIncludePath(ceres.path("internal/ceres/miniglog"));
    wasm_lib_mod.addIncludePath(eigen.path("."));
    wasm_lib_mod.addIncludePath(b.path("include"));

    const wasm_mod = b.createModule(.{
        .root_source_file = b.path("src/wasm.zig"),
        .target = wasm_target,
        .optimize = optimize,
    });

    wasm_mod.addImport("ozric_lib", wasm_lib_mod);

    const wasm_exe = b.addExecutable(.{
        .name = "ozric_wasm",
        .root_module = wasm_mod,
    });

    // WASM-specific flags (more restrictive than native)
    const wasm_ceres_flags = &.{
        "-std=c++17",
        "-DCERES_NO_SUITESPARSE",
        "-DCERES_NO_CXSPARSE",
        "-DCERES_EXPORT=",
        "-DCERES_NO_EXPORT=",
        "-DCERES_NO_PROTOCOL_BUFFERS",
        "-DCERES_NO_THREADS",
        "-DCERES_NO_CUDA",
        "-DCERES_NO_LAPACK",
        "-DCERES_METIS_VERSION=\"5.1.0\"",
        "-DMINIGLOG",
        "-DCERES_RESTRICT_SCHUR_SPECIALIZATION",
        "-DCERES_NO_ACCELERATE_SPARSE",
        "-DCERES_NO_CUSTOM_BLAS",
        "-fno-exceptions", // Disable C++ exceptions for WASM
        "-fno-rtti", // Disable runtime type info for smaller binary
        "-D__wasm__", // Define __wasm__ for our threading stub detection
    };

    // Ultra-minimal WASM sources (just test version number for now)
    const wasm_ceres_sources = [_][]const u8{
        "stringprintf.cc",
        "wall_time.cc",
        "types.cc",
        "cost_function.cc",
        "loss_function.cc",
        "is_close.cc",
    };

    // Add C++ sources directly to WASM executable
    wasm_exe.addCSourceFiles(.{
        .root = b.path("src"),
        .files = &.{"solver.cc"},
        .flags = wasm_ceres_flags,
    });

    wasm_exe.addCSourceFiles(.{
        .root = ceres.path("internal/ceres"),
        .files = &wasm_ceres_sources,
        .flags = wasm_ceres_flags,
    });

    wasm_exe.addCSourceFiles(.{
        .root = ceres.path("internal/ceres/miniglog/glog"),
        .files = &.{"logging.cc"},
        .flags = wasm_ceres_flags,
    });

    // Add include paths to WASM executable
    wasm_exe.addIncludePath(ceres.path("include"));
    wasm_exe.addIncludePath(ceres.path("."));
    wasm_exe.addIncludePath(ceres.path("internal"));
    wasm_exe.addIncludePath(ceres.path("config"));
    wasm_exe.addIncludePath(ceres.path("internal/ceres/miniglog"));
    wasm_exe.addIncludePath(eigen.path("."));
    wasm_exe.addIncludePath(b.path("include"));

    wasm_exe.linkLibC();
    wasm_exe.linkLibCpp();

    b.installArtifact(wasm_exe);

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
}
