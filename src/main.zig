const std = @import("std");
const ozric = @import("ozric");

/// Helper function to create file paths in the output directory
fn buildOutputPath(allocator: std.mem.Allocator, output_dir: []const u8, filename: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, "{s}/{s}", .{ output_dir, filename });
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.debug.print("🧪 Ozric - Ornstein-Zernike Equation Solver", .{});
    std.debug.print("", .{});

    // Parse command line arguments
    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    // Parse arguments: ozric [example] [--output-dir <dir>]
    var example: []const u8 = "hard-sphere";
    var output_dir: []const u8 = "out";

    var i: usize = 1;
    while (i < args.len) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "help") or std.mem.eql(u8, arg, "--help")) {
            printUsage();
            return;
        } else if (std.mem.eql(u8, arg, "--output-dir") or std.mem.eql(u8, arg, "-o")) {
            if (i + 1 >= args.len) {
                std.log.err("Error: {s} requires a directory argument", .{arg});
                return;
            }
            i += 1;
            output_dir = args[i];
        } else {
            // First non-flag argument is the example
            example = arg;
        }
        i += 1;
    }

    // Run the selected example
    if (std.mem.eql(u8, example, "hard-sphere")) {
        try runHardSphereExample(allocator, output_dir);
    } else if (std.mem.eql(u8, example, "lennard-jones")) {
        try runLennardJonesExample(allocator, output_dir);
    } else if (std.mem.eql(u8, example, "comparison")) {
        try runMethodComparison(allocator, output_dir);
    } else {
        std.log.err("Unknown example: '{s}'. Use 'help' to see available examples.", .{example});
        return;
    }
}

fn printUsage() void {
    std.debug.print("Usage: ozric [example] [options]\n", .{});
    std.debug.print("\n", .{});
    std.debug.print("Available examples:\n", .{});
    std.debug.print("  hard-sphere    - Hard sphere potential with PY closure (default)\n", .{});
    std.debug.print("  lennard-jones  - Lennard-Jones potential with PY closure\n", .{});
    std.debug.print("  comparison     - Compare different solver methods\n", .{});
    std.debug.print("  help           - Show this help message\n", .{});
    std.debug.print("\n", .{});
    std.debug.print("Options:\n", .{});
    std.debug.print("  -o, --output-dir <dir>  Output directory (default: 'out')\n", .{});
    std.debug.print("\n", .{});
    std.debug.print("Examples:\n", .{});
    std.debug.print("  ozric hard-sphere\n", .{});
    std.debug.print("  ozric lennard-jones --output-dir results\n", .{});
    std.debug.print("  ozric comparison -o /tmp/data\n", .{});
}

fn runHardSphereExample(allocator: std.mem.Allocator, output_dir: []const u8) !void {
    std.log.info("🔴 Hard Sphere System Example", .{});
    std.log.info("", .{});

    const grid = ozric.GridParams.init(256, 8.0);
    var solver = try ozric.Solver.init(allocator, grid);
    defer solver.deinit();

    // Hard sphere parameters: moderate density
    const density = 0.7;
    const temperature = 1.0;
    const sigma = 1.0;

    std.log.info("System parameters:", .{});
    std.log.info("  Density (ρ):     {d:.3}", .{density});
    std.log.info("  Temperature (T): {d:.3}", .{temperature});
    std.log.info("  Diameter (σ):    {d:.3}", .{sigma});
    std.log.info("", .{});

    solver.initHardSphere(density, temperature, sigma);

    // Solve with conservative parameters
    const params = ozric.ConvergenceParams{
        .method = .adaptive_mixing,
        .initial_mixing = 0.1,
        .min_mixing = 0.01,
        .max_mixing = 0.3,
    };
    solver.setConvergenceParams(params);

    std.log.info("🔧 Solving OZ equation...", .{});
    try solver.solve(.percus_yevick, .simple_convolution, 100, 1e-5);

    // Display results
    ozric.export_data.displayRadialFunctions(&solver, 15);

    // ASCII plot
    try ozric.export_data.plotRadialDistribution(&solver, 60, 15);

    // Export data
    std.log.info("💾 Exporting data...", .{});

    // Ensure output directory exists
    std.fs.cwd().makeDir(output_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };
    const csv_path = try buildOutputPath(allocator, output_dir, "hard_sphere_rdf.csv");
    defer allocator.free(csv_path);
    try ozric.export_data.exportRadialData(allocator, &solver, csv_path, csv_config);

    const json_config = ozric.export_data.ExportConfig{ .format = .json, .precision = 6 };
    const json_path = try buildOutputPath(allocator, output_dir, "hard_sphere_rdf.json");
    defer allocator.free(json_path);
    try ozric.export_data.exportRadialData(allocator, &solver, json_path, json_config);

    const summary_path = try buildOutputPath(allocator, output_dir, "hard_sphere_summary.txt");
    defer allocator.free(summary_path);
    try ozric.export_data.exportSummaryStats(&solver, summary_path);

    std.log.info("✅ Files saved:", .{});
    std.log.info("  - {s}/hard_sphere_rdf.csv     (CSV data)", .{output_dir});
    std.log.info("  - {s}/hard_sphere_rdf.json    (JSON data)", .{output_dir});
    std.log.info("  - {s}/hard_sphere_summary.txt (Statistics)", .{output_dir});
}

fn runLennardJonesExample(allocator: std.mem.Allocator, output_dir: []const u8) !void {
    std.log.info("⚛️  Lennard-Jones System Example", .{});
    std.log.info("", .{});

    const grid = ozric.GridParams.init(256, 8.0);
    var solver = try ozric.Solver.init(allocator, grid);
    defer solver.deinit();

    // LJ parameters: liquid-like conditions
    const density = 0.8;
    const temperature = 1.5;
    const epsilon = 1.0;
    const sigma = 1.0;

    std.log.info("System parameters:", .{});
    std.log.info("  Density (ρ):     {d:.3}", .{density});
    std.log.info("  Temperature (T): {d:.3}", .{temperature});
    std.log.info("  Energy (ε):      {d:.3}", .{epsilon});
    std.log.info("  Length (σ):      {d:.3}", .{sigma});
    std.log.info("", .{});

    solver.initLennardJones(density, temperature, epsilon, sigma);

    // Use conservative convergence for stability
    const params = ozric.ConvergenceParams{
        .method = .adaptive_mixing,
        .initial_mixing = 0.05,
        .min_mixing = 0.001,
        .max_mixing = 0.2,
        .error_increase_factor = 1.2,
    };
    solver.setConvergenceParams(params);

    std.log.info("🔧 Solving OZ equation...", .{});
    try solver.solve(.percus_yevick, .simple_convolution, 150, 1e-4);

    // Display results
    ozric.export_data.displayRadialFunctions(&solver, 15);

    // ASCII plot
    try ozric.export_data.plotRadialDistribution(&solver, 60, 15);

    // Export data
    std.log.info("💾 Exporting data...", .{});

    // Ensure output directory exists
    std.fs.cwd().makeDir(output_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };
    const csv_path = try buildOutputPath(allocator, output_dir, "lennard_jones_rdf.csv");
    defer allocator.free(csv_path);
    try ozric.export_data.exportRadialData(allocator, &solver, csv_path, csv_config);

    const json_config = ozric.export_data.ExportConfig{ .format = .json, .precision = 6 };
    const json_path = try buildOutputPath(allocator, output_dir, "lennard_jones_rdf.json");
    defer allocator.free(json_path);
    try ozric.export_data.exportRadialData(allocator, &solver, json_path, json_config);

    const summary_path = try buildOutputPath(allocator, output_dir, "lennard_jones_summary.txt");
    defer allocator.free(summary_path);
    try ozric.export_data.exportSummaryStats(&solver, summary_path);

    std.log.info("✅ Files saved:", .{});
    std.log.info("  - {s}/lennard_jones_rdf.csv     (CSV data)", .{output_dir});
    std.log.info("  - {s}/lennard_jones_rdf.json    (JSON data)", .{output_dir});
    std.log.info("  - {s}/lennard_jones_summary.txt (Statistics)", .{output_dir});
}

fn runMethodComparison(allocator: std.mem.Allocator, output_dir: []const u8) !void {
    std.log.info("⚖️  Solver Method Comparison", .{});
    std.log.info("", .{});

    const grid = ozric.GridParams.init(128, 6.0);

    // Test system: hard spheres at moderate density
    const density = 0.6;
    const temperature = 1.0;
    const sigma = 1.0;

    std.log.info("Comparing simple convolution vs FFT-based solver", .{});
    std.log.info("System: Hard spheres, ρ={d:.2}, T={d:.2}", .{ density, temperature });
    std.log.info("", .{});

    // Simple convolution method
    std.log.info("🔄 Testing simple convolution method...", .{});
    var solver_simple = try ozric.Solver.init(allocator, grid);
    defer solver_simple.deinit();

    solver_simple.initHardSphere(density, temperature, sigma);
    const start_simple = std.time.milliTimestamp();
    try solver_simple.solve(.percus_yevick, .simple_convolution, 50, 1e-4);
    const time_simple = std.time.milliTimestamp() - start_simple;

    // FFT-based method
    std.log.info("🌊 Testing FFT-based method...", .{});
    var solver_fft = try ozric.Solver.init(allocator, grid);
    defer solver_fft.deinit();

    solver_fft.initHardSphere(density, temperature, sigma);
    const start_fft = std.time.milliTimestamp();
    try solver_fft.solve(.percus_yevick, .fft_based, 50, 1e-4);
    const time_fft = std.time.milliTimestamp() - start_fft;

    // Compare results
    std.log.info("⏱️  Performance comparison:", .{});
    std.log.info("  Simple convolution: {}ms", .{time_simple});
    std.log.info("  FFT-based:          {}ms", .{time_fft});
    std.log.info("", .{});

    // Show both results briefly
    std.log.info("📊 Simple convolution g(r):", .{});
    ozric.export_data.displayRadialFunctions(&solver_simple, 8);

    std.log.info("📊 FFT-based g(r):", .{});
    ozric.export_data.displayRadialFunctions(&solver_fft, 8);

    // Export comparison data
    std.fs.cwd().makeDir(output_dir) catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };

    const simple_path = try buildOutputPath(allocator, output_dir, "comparison_simple.csv");
    defer allocator.free(simple_path);
    try ozric.export_data.exportRadialData(allocator, &solver_simple, simple_path, csv_config);

    const fft_path = try buildOutputPath(allocator, output_dir, "comparison_fft.csv");
    defer allocator.free(fft_path);
    try ozric.export_data.exportRadialData(allocator, &solver_fft, fft_path, csv_config);

    std.log.info("✅ Comparison files saved:", .{});
    std.log.info("  - {s}/comparison_simple.csv (Simple convolution)", .{output_dir});
    std.log.info("  - {s}/comparison_fft.csv    (FFT-based)", .{output_dir});
}
