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
    lib.linkSystemLibrary("nlopt");

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
    exe.linkSystemLibrary("nlopt");

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
    lib_unit_tests.linkSystemLibrary("nlopt");

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
    exe_unit_tests.linkSystemLibrary("nlopt");

    const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);
    test_step.dependOn(&run_exe_unit_tests.step);

    // WebAssembly target
    const wasm_target = b.resolveTargetQuery(.{
        .cpu_arch = .wasm32,
        .os_tag = .freestanding,
    });

    const wasm_lib = b.addStaticLibrary(.{
        .name = "ozric",
        .root_source_file = b.path("src/wasm.zig"),
        .target = wasm_target,
        .optimize = optimize,
    });

    // Also create a WebAssembly module (executable)
    const wasm_exe = b.addExecutable(.{
        .name = "ozric",
        .root_source_file = b.path("src/wasm.zig"),
        .target = wasm_target,
        .optimize = optimize,
    });

    const wasm_step = b.step("wasm", "Build WebAssembly library and module");
    wasm_step.dependOn(&b.addInstallArtifact(wasm_lib, .{}).step);
    wasm_step.dependOn(&b.addInstallArtifact(wasm_exe, .{}).step);
}
