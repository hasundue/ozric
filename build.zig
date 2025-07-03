const std = @import("std");
const ceres = @import("src/build/ceres.zig");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Dependencies
    const ceres_dep = b.dependency("ceres", .{
        .target = target,
        .optimize = optimize,
    });
    const eigen = b.dependency("eigen", .{
        .target = target,
        .optimize = optimize,
    });

    // Library
    const lib = blk: {
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

        // Add Ceres-solver support using the extracted module
        ceres.addCeresSupport(b, lib, ceres_dep, eigen);

        b.installArtifact(lib);
        break :blk .{ .artifact = lib, .module = lib_mod };
    };

    // Executable
    const exe = blk: {
        const exe_mod = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        });

        exe_mod.addImport("ozric_lib", lib.module);

        const exe = b.addExecutable(.{
            .name = "ozric",
            .root_module = exe_mod,
        });

        b.installArtifact(exe);
        break :blk .{ .artifact = exe, .module = exe_mod };
    };

    // Run step
    {
        const run_cmd = b.addRunArtifact(exe.artifact);
        run_cmd.step.dependOn(b.getInstallStep());

        if (b.args) |args| {
            run_cmd.addArgs(args);
        }

        const run_step = b.step("run", "Run the app");
        run_step.dependOn(&run_cmd.step);
    }

    // Test step
    {
        const lib_unit_tests = b.addTest(.{
            .root_module = lib.module,
        });
        const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

        const exe_unit_tests = b.addTest(.{
            .root_module = exe.module,
        });
        const run_exe_unit_tests = b.addRunArtifact(exe_unit_tests);

        const test_step = b.step("test", "Run unit tests");
        test_step.dependOn(&run_lib_unit_tests.step);
        test_step.dependOn(&run_exe_unit_tests.step);
    }

    // WASM step
    {
        const wasm_step = b.step("wasm", "Build WASM using Emscripten");

        // Get all source files for WASM compilation
        var all_sources = ceres.getCeresSources(ceres_dep, b.allocator) catch @panic("OOM");

        // Add our project-specific C++ wrapper
        all_sources.append(.{
            .path = "src/solver.cc",
            .base_dir = b.path("."),
        }) catch @panic("OOM");

        // Compile source files to object files
        var object_files = std.ArrayList([]const u8).init(b.allocator);
        var compile_steps = std.ArrayList(*std.Build.Step).init(b.allocator);

        for (all_sources.items) |source_info| {
            const compile_cmd = b.addSystemCommand(&ceres.emcc_compile_command);

            // Set environment variables
            compile_cmd.setEnvironmentVariable("EMCC_DEBUG", "1");
            compile_cmd.setEnvironmentVariable("CCACHE_DIR", ".zig-cache/ccache");
            compile_cmd.setEnvironmentVariable("CCACHE_IGNOREOPTIONS", "-v -g* -W* --verbose");

            // Add include paths using the ceres module
            ceres.addEmccIncludePaths(compile_cmd, b, ceres_dep, eigen);

            // Input source file
            const source_path = b.fmt("{s}/{s}", .{ source_info.base_dir.getPath(b), source_info.path });
            compile_cmd.addArg(source_path);

            // Output object file with directory structure to avoid collisions
            const source_basename = std.fs.path.basename(source_info.path);
            const base_name = source_basename[0 .. source_basename.len - 3]; // Remove .cc/.cpp extension

            // Organize object files by source: ours in obj/src/, ceres in obj/ceres/
            const object_path = if (std.mem.startsWith(u8, source_info.path, "src/"))
                b.fmt("zig-out/obj/src/{s}.o", .{base_name})
            else
                b.fmt("zig-out/obj/ceres/{s}/{s}.o", .{ std.fs.path.dirname(source_info.path) orelse "", base_name });
            compile_cmd.addArg("-o");
            compile_cmd.addArg(object_path);

            object_files.append(object_path) catch @panic("OOM");
            compile_steps.append(&compile_cmd.step) catch @panic("OOM");
        }

        // Final linking step
        const link_cmd = b.addSystemCommand(&.{
            "emcc",
            "-pthread",
            "-sASYNCIFY",
            "-sPTHREAD_POOL_SIZE=4",
            "-sINITIAL_MEMORY=167772160", // 160MB
            "-sUSE_OFFSET_CONVERTER",
            "-sEXPORTED_FUNCTIONS=['_test_ceres','_run_hello_world']",
            "-sEXPORTED_RUNTIME_METHODS=ccall,cwrap,PThread",
            "-O0",
            "--cache",
            ".zig-cache/emcc",
        });

        link_cmd.setEnvironmentVariable("EMCC_DEBUG", "1");

        // Add all object files
        for (object_files.items) |obj_file| {
            link_cmd.addArg(obj_file);
        }

        // Output file
        link_cmd.addArg("-o");
        link_cmd.addArg("zig-out/bin/ozric_wasm.js");

        // Make link_cmd depend on all compile steps
        for (compile_steps.items) |compile_step| {
            link_cmd.step.dependOn(compile_step);
        }

        wasm_step.dependOn(&link_cmd.step);
    }

    // Compile commands step
    {
        const compile_commands_step = b.step("compile-commands", "Generate compile_commands.json for clangd");

        const gen_compile_commands = b.addWriteFiles();

        // Build the compile command using the same flags and includes as the actual build
        var compile_command = std.ArrayList(u8).init(b.allocator);
        const writer = compile_command.writer();

        // Start with clang++
        writer.print("clang++", .{}) catch @panic("OOM");

        // Add all ceres flags
        for (ceres.ceres_flags) |flag| {
            writer.print(" {s}", .{flag}) catch @panic("OOM");
        }

        // Add include paths using the same function as addCeresSupport
        ceres.addIncludePathsToWriter(writer, b, ceres_dep, eigen) catch @panic("OOM");
        writer.print(" -c src/solver.cc", .{}) catch @panic("OOM");

        // Generate the JSON content with proper escaping
        // Need to escape quotes in the command string for valid JSON
        var escaped_command = std.ArrayList(u8).init(b.allocator);
        const escape_writer = escaped_command.writer();
        for (compile_command.items) |char| {
            if (char == '"') {
                escape_writer.print("\\\"", .{}) catch @panic("OOM");
            } else {
                escape_writer.print("{c}", .{char}) catch @panic("OOM");
            }
        }

        const json_content = b.fmt(
            \\[
            \\  {{
            \\    "directory": "{s}",
            \\    "command": "{s}",
            \\    "file": "src/solver.cc"
            \\  }}
            \\]
            \\
        , .{ b.path("").getPath(b), escaped_command.items });

        const compile_commands_file = gen_compile_commands.add("compile_commands.json", json_content);

        const install_compile_commands = b.addInstallFile(compile_commands_file, "../compile_commands.json");

        compile_commands_step.dependOn(&install_compile_commands.step);
    }
}
