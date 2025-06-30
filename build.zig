const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const ceres = b.dependency("ceres", .{
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
        .name = "default",
        .root_module = lib_mod,
    });

    lib.addCSourceFiles(.{
        .root = b.path("src"),
        .files = &.{"solver.cc"},
    });

    // Add all ceres source files for static compilation
    const ceres_sources = [_][]const u8{
        "sparse_cholesky.cc",              "file.cc",                         "problem_impl.cc",                     "cgnr_solver.cc",
        "levenberg_marquardt_strategy.cc", "dense_normal_cholesky_solver.cc", "block_random_access_dense_matrix.cc", "normal_prior.cc",
        "solver.cc",                       "problem.cc",                      "covariance_impl.cc",
        "gradient_problem.cc",
        // Add more core files as needed - this is a minimal set
    };

    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres"),
        .files = &ceres_sources,
        .flags = &.{ "-std=c++17", "-DCERES_NO_SUITESPARSE", "-DCERES_NO_CXSPARSE" },
    });

    // Add ceres include directories
    lib.addIncludePath(ceres.path("include"));
    lib.addIncludePath(ceres.path(".")); // For internal headers
    lib.linkLibC();
    lib.linkLibCpp();

    lib.installHeadersDirectory(ceres.path("include"), "", .{
        .include_extensions = &.{".h"},
    });

    b.installArtifact(lib);

    const exe_mod = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    exe_mod.addImport("default_lib", lib_mod);

    const exe = b.addExecutable(.{
        .name = "default",
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
