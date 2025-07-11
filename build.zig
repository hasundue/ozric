const std = @import("std");
const b_ceres = @import("src/build/ceres.zig");

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

    // Zig bindings to Ceres
    const ceres = blk: {
        const mod = b.createModule(.{
            .root_source_file = b.path("src/ceres.zig"),
            .target = target,
            .optimize = optimize,
        });

        const lib = b.addLibrary(.{
            .linkage = .static,
            .name = "ceres",
            .root_module = mod,
        });

        lib.addCSourceFiles(.{
            .root = b.path("src"),
            .files = &.{"ceres.cc"},
        });

        b_ceres.addCeresSupport(b, lib, ceres_dep, eigen);

        const tests = b.addTest(.{
            .root_module = mod,
        });

        break :blk .{
            .mod = mod,
            .lib = lib,
            .tests = tests,
        };
    };

    // Build and test Ceres bindings
    {
        const step = b.step("ceres", "Test Ceres bindings");
        const run_tests = b.addRunArtifact(ceres.tests);
        step.dependOn(&run_tests.step);
    }

    // Library
    const lib = blk: {
        const mod = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
            .optimize = optimize,
        });

        const lib = b.addLibrary(.{
            .linkage = .static,
            .name = "ozric",
            .root_module = mod,
        });

        const tests = b.addTest(.{
            .root_module = mod,
        });

        break :blk .{ .artifact = lib, .module = mod, .tests = tests };
    };

    // Test Zig library
    {
        const step = b.step("test", "Run unit tests");
        const run = b.addRunArtifact(lib.tests);
        step.dependOn(&run.step);
    }

    // Integration tests
    {
        const step = b.step("test-integration", "Run integration tests");
        const integration_tests = b.addTest(.{
            .root_source_file = b.path("tests/convolution.zig"),
            .target = target,
            .optimize = optimize,
        });
        integration_tests.root_module.addImport("convolution", lib.module);
        const run = b.addRunArtifact(integration_tests);
        step.dependOn(&run.step);
    }

    // WASM step
    {
        const wasm_step = b.step("wasm", "Build WASM using Emscripten");

        // Get all source files for WASM compilation
        var all_sources = b_ceres.getCeresSources(ceres_dep, b.allocator) catch @panic("OOM");

        // Add our project-specific C++ wrapper
        all_sources.append(.{
            .path = "src/solver.cc",
            .base_dir = b.path("."),
        }) catch @panic("OOM");

        // Compile source files to object files
        var object_files = std.ArrayList([]const u8).init(b.allocator);
        var compile_steps = std.ArrayList(*std.Build.Step).init(b.allocator);

        for (all_sources.items) |source_info| {
            const compile_cmd = b.addSystemCommand(&b_ceres.emcc_compile_command);

            // Set environment variables
            compile_cmd.setEnvironmentVariable("EMCC_DEBUG", "1");
            compile_cmd.setEnvironmentVariable("CCACHE_DIR", ".zig-cache/ccache");
            compile_cmd.setEnvironmentVariable("CCACHE_IGNOREOPTIONS", "-v -g* -W* --verbose");

            // Add include paths using the ceres module
            b_ceres.addEmccIncludePaths(compile_cmd, b, ceres_dep, eigen);

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
            "-sEXPORTED_FUNCTIONS=['_test_ceres','_solve']",
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
        for (b_ceres.ceres_flags) |flag| {
            writer.print(" {s}", .{flag}) catch @panic("OOM");
        }

        // Add include paths using the same function as addCeresSupport
        b_ceres.addIncludePathsToWriter(writer, b, ceres_dep, eigen) catch @panic("OOM");
        writer.print(" -c src/solver.cc", .{}) catch @panic("OOM");

        // Serialize to JSON using std.json
        var json_content = std.ArrayList(u8).init(b.allocator);
        const json_writer = json_content.writer();
        std.json.stringify(.{
            .{
                .directory = b.path("").getPath(b),
                .command = compile_command.items,
                .file = "src/solver.cc",
            },
        }, .{}, json_writer) catch @panic("JSON serialization failed");

        const compile_commands_file = gen_compile_commands.add("compile_commands.json", json_content.items);
        const install_compile_commands = b.addInstallFile(compile_commands_file, "../compile_commands.json");
        compile_commands_step.dependOn(&install_compile_commands.step);
    }

    // Fish completions step
    {
        const fish_completions_step = b.step("fish-completions", "Generate fish shell completions");

        const gen_fish_completions = b.addWriteFiles();

        // Generate fish completion content dynamically from build steps
        var fish_content = std.ArrayList(u8).init(b.allocator);
        const writer = fish_content.writer();

        // Header
        writer.print(
            \\# Ozric project-specific Zig build completions
            \\# Generated automatically by 'zig build fish-completions'
            \\# 
            \\# To use: source completions/zig.fish
            \\
            \\# Add our custom build steps (supplements built-in fish zig completions)
            \\
        , .{}) catch @panic("OOM");

        // Add completions for all registered build steps
        var step_iterator = b.top_level_steps.iterator();
        while (step_iterator.next()) |entry| {
            const step_name = entry.key_ptr.*;
            const step_info = entry.value_ptr.*;

            // Skip standard steps that fish already knows about
            if (std.mem.eql(u8, step_name, "install") or
                std.mem.eql(u8, step_name, "uninstall") or
                std.mem.eql(u8, step_name, "run") or
                std.mem.eql(u8, step_name, "test"))
            {
                continue;
            }

            // Add completion for this step
            writer.print("complete -c zig -n '__fish_seen_subcommand_from build' -f -a '{s}' -d '{s}'\n", .{ step_name, step_info.description }) catch @panic("OOM");
        }

        const fish_file = gen_fish_completions.add("zig.fish", fish_content.items);

        const install_fish_completions = b.addInstallFile(fish_file, "../completions/zig.fish");

        fish_completions_step.dependOn(&install_fish_completions.step);
    }
}
