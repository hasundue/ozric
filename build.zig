const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Library
    const lib = b.addStaticLibrary(.{
        .name = "ozric",
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Link C math libraries
    lib.linkLibC();
    lib.linkSystemLibrary("blas");
    lib.linkSystemLibrary("lapack");
    lib.linkSystemLibrary("fftw3");

    b.installArtifact(lib);

    // Executable
    const exe = b.addExecutable(.{
        .name = "ozric",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    exe.root_module.addImport("ozric", &lib.root_module.*);
    exe.linkLibC();
    exe.linkSystemLibrary("blas");
    exe.linkSystemLibrary("lapack");
    exe.linkSystemLibrary("fftw3");

    b.installArtifact(exe);

    // Run step
    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    // Tests
    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    lib_unit_tests.linkLibC();
    lib_unit_tests.linkSystemLibrary("blas");
    lib_unit_tests.linkSystemLibrary("lapack");
    lib_unit_tests.linkSystemLibrary("fftw3");

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    const exe_unit_tests = b.addTest(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });

    exe_unit_tests.root_module.addImport("ozric", &lib.root_module.*);
    exe_unit_tests.linkLibC();
    exe_unit_tests.linkSystemLibrary("blas");
    exe_unit_tests.linkSystemLibrary("lapack");
    exe_unit_tests.linkSystemLibrary("fftw3");

    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
    test_step.dependOn(&run_exe_unit_tests.step);
}