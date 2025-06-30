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

        // Dense matrix support
        "dense_sparse_matrix.cc",

        // Parallel utilities
        "parallel_utils.cc",

        // Context
        "context.cc",
        "context_impl.cc",
    };

    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres"),
        .files = &essential_ceres_sources,
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

    lib.linkLibC();
    lib.linkLibCpp();

    b.installArtifact(lib);

    // TODO: Re-enable WASM when ceres compilation is stable
    // const wasm_mod = b.createModule(.{
    //     .root_source_file = b.path("src/wasm.zig"),
    //     .target = b.resolveTargetQuery(.{
    //         .cpu_arch = .wasm32,
    //         .os_tag = .wasi,
    //     }),
    //     .optimize = optimize,
    // });

    // wasm_mod.addImport("ozric_lib", lib_mod);

    // const wasm_lib = b.addLibrary(.{
    //     .linkage = .static,
    //     .name = "libozric_wasm",
    //     .root_module = wasm_mod,
    // });

    // b.installArtifact(wasm_lib);

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
