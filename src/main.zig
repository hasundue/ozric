const std = @import("std");
const ozric = @import("ozric");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    std.log.info("🧪 Ozric - Ornstein-Zernike Equation Solver", .{});
    std.log.info("", .{});

    // Parse command line arguments for different examples
    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    const example = if (args.len > 1) args[1] else "hard-sphere";

    if (std.mem.eql(u8, example, "help") or std.mem.eql(u8, example, "--help")) {
        printUsage();
        return;
    }

    // Run the selected example
    if (std.mem.eql(u8, example, "hard-sphere")) {
        try runHardSphereExample(allocator);
    } else if (std.mem.eql(u8, example, "lennard-jones")) {
        try runLennardJonesExample(allocator);
    } else if (std.mem.eql(u8, example, "comparison")) {
        try runMethodComparison(allocator);
    } else {
        std.log.err("Unknown example: '{s}'. Use 'help' to see available examples.", .{example});
        return;
    }
}

fn printUsage() void {
    std.log.info("Usage: ozric [example]", .{});
    std.log.info("", .{});
    std.log.info("Available examples:", .{});
    std.log.info("  hard-sphere    - Hard sphere potential with PY closure (default)", .{});
    std.log.info("  lennard-jones  - Lennard-Jones potential with PY closure", .{});
    std.log.info("  comparison     - Compare different solver methods", .{});
    std.log.info("  help           - Show this help message", .{});
    std.log.info("", .{});
    std.log.info("Output files will be saved to the 'out/' directory.", .{});
}

fn runHardSphereExample(allocator: std.mem.Allocator) !void {
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
    std.fs.cwd().makeDir("out") catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };
    try ozric.export_data.exportRadialData(allocator, &solver, "out/hard_sphere_rdf.csv", csv_config);

    const json_config = ozric.export_data.ExportConfig{ .format = .json, .precision = 6 };
    try ozric.export_data.exportRadialData(allocator, &solver, "out/hard_sphere_rdf.json", json_config);

    try ozric.export_data.exportSummaryStats(&solver, "out/hard_sphere_summary.txt");

    std.log.info("✅ Files saved:", .{});
    std.log.info("  - out/hard_sphere_rdf.csv     (CSV data)", .{});
    std.log.info("  - out/hard_sphere_rdf.json    (JSON data)", .{});
    std.log.info("  - out/hard_sphere_summary.txt (Statistics)", .{});
}

fn runLennardJonesExample(allocator: std.mem.Allocator) !void {
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
    std.fs.cwd().makeDir("out") catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };
    try ozric.export_data.exportRadialData(allocator, &solver, "out/lennard_jones_rdf.csv", csv_config);

    const json_config = ozric.export_data.ExportConfig{ .format = .json, .precision = 6 };
    try ozric.export_data.exportRadialData(allocator, &solver, "out/lennard_jones_rdf.json", json_config);

    try ozric.export_data.exportSummaryStats(&solver, "out/lennard_jones_summary.txt");

    std.log.info("✅ Files saved:", .{});
    std.log.info("  - out/lennard_jones_rdf.csv     (CSV data)", .{});
    std.log.info("  - out/lennard_jones_rdf.json    (JSON data)", .{});
    std.log.info("  - out/lennard_jones_summary.txt (Statistics)", .{});
}

fn runMethodComparison(allocator: std.mem.Allocator) !void {
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
    std.fs.cwd().makeDir("out") catch |err| switch (err) {
        error.PathAlreadyExists => {},
        else => return err,
    };

    const csv_config = ozric.export_data.ExportConfig{ .format = .csv, .precision = 6 };
    try ozric.export_data.exportRadialData(allocator, &solver_simple, "out/comparison_simple.csv", csv_config);
    try ozric.export_data.exportRadialData(allocator, &solver_fft, "out/comparison_fft.csv", csv_config);

    std.log.info("✅ Comparison files saved:", .{});
    std.log.info("  - out/comparison_simple.csv (Simple convolution)", .{});
    std.log.info("  - out/comparison_fft.csv    (FFT-based)", .{});
}
