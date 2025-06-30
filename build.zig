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
        .name = "default",
        .root_module = lib_mod,
    });

    lib.addCSourceFiles(.{
        .root = b.path("src"),
        .files = &.{"solver.cc"},
    });

    // Add minimal ceres source files with MINIGLOG for 2.2.0
    const ceres_sources = [_][]const u8{
        "file.cc", // Simple utility file that now uses glog instead of absl
    };

    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres"),
        .files = &ceres_sources,
        .flags = &.{
            "-std=c++17",
            "-DCERES_NO_SUITESPARSE",
            "-DCERES_NO_CXSPARSE",
            "-DCERES_EXPORT=", // Define CERES_EXPORT as empty
            "-DCERES_NO_EXPORT=", // Define CERES_NO_EXPORT as empty
            "-DCERES_NO_PROTOCOL_BUFFERS",
            "-DCERES_NO_THREADS",
            "-DMINIGLOG", // Use minimal logging instead of full glog
        },
    });

    // Add miniglog source file
    lib.addCSourceFiles(.{
        .root = ceres.path("internal/ceres/miniglog/glog"),
        .files = &.{"logging.cc"},
        .flags = &.{
            "-std=c++17",
            "-DCERES_EXPORT=", // Define CERES_EXPORT as empty for miniglog too
            "-DCERES_NO_EXPORT=", // Define CERES_NO_EXPORT as empty for miniglog too
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
