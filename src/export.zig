const std = @import("std");
const grid = @import("grid.zig");
const solver = @import("solver.zig");

const RadialFunction = grid.RadialFunction;
const Solver = solver.Solver;

/// Export format options
pub const ExportFormat = enum {
    csv,
    json,
    tsv, // Tab-separated values
};

/// Export configuration
pub const ExportConfig = struct {
    format: ExportFormat = .csv,
    include_header: bool = true,
    precision: u8 = 6, // Number of decimal places
    delimiter: []const u8 = ",", // For CSV/TSV
};

/// Export error types
pub const ExportError = error{
    FileCreationFailed,
    WriteError,
    InvalidFormat,
};

/// Export radial distribution data to file
pub fn exportRadialData(
    allocator: std.mem.Allocator,
    oz_solver: *const Solver,
    file_path: []const u8,
    config: ExportConfig,
) !void {
    const file = std.fs.cwd().createFile(file_path, .{}) catch |err| switch (err) {
        error.PathAlreadyExists => try std.fs.cwd().openFile(file_path, .{ .mode = .write_only }),
        else => return ExportError.FileCreationFailed,
    };
    defer file.close();

    const writer = file.writer();

    switch (config.format) {
        .csv => try writeCSV(writer, oz_solver, config),
        .json => try writeJSON(allocator, writer, oz_solver, config),
        .tsv => {
            var tsv_config = config;
            tsv_config.delimiter = "\t";
            try writeCSV(writer, oz_solver, tsv_config);
        },
    }
}

/// Write data in CSV format
fn writeCSV(
    writer: anytype,
    oz_solver: *const Solver,
    config: ExportConfig,
) !void {
    // Write header if requested
    if (config.include_header) {
        try writer.print("r{s}g_r{s}h_r{s}c_r\n", .{ config.delimiter, config.delimiter, config.delimiter });
    }

    // Write data rows
    for (0..oz_solver.grid.n_points) |i| {
        const r = oz_solver.g_r.getRadius(i);
        const g_r = oz_solver.g_r.values[i];
        const h_r = oz_solver.h_r.values[i];
        const c_r = oz_solver.c_r.values[i];

        try writer.print("{d:.6}{s}{d:.6}{s}{d:.6}{s}{d:.6}\n", .{
            r,   config.delimiter,
            g_r, config.delimiter,
            h_r, config.delimiter,
            c_r,
        });
    }
}

/// Write data in JSON format
fn writeJSON(
    _: std.mem.Allocator,
    writer: anytype,
    oz_solver: *const Solver,
    _: ExportConfig,
) !void {
    // Build JSON structure
    try writer.writeAll("{\n");
    try writer.writeAll("  \"system_info\": {\n");
    try writer.print("    \"density\": {d:.6},\n", .{oz_solver.density});
    try writer.print("    \"temperature\": {d:.6},\n", .{oz_solver.temperature});
    try writer.print("    \"beta\": {d:.6},\n", .{oz_solver.beta});
    try writer.print("    \"n_points\": {},\n", .{oz_solver.grid.n_points});
    try writer.print("    \"r_max\": {d:.6},\n", .{oz_solver.grid.r_max});
    try writer.print("    \"dr\": {d:.6}\n", .{oz_solver.grid.dr});
    try writer.writeAll("  },\n");

    // Write correlation function data
    try writer.writeAll("  \"data\": [\n");
    for (0..oz_solver.grid.n_points) |i| {
        const r = oz_solver.g_r.getRadius(i);
        const g_r = oz_solver.g_r.values[i];
        const h_r = oz_solver.h_r.values[i];
        const c_r = oz_solver.c_r.values[i];

        try writer.writeAll("    {\n");
        try writer.print("      \"r\": {d:.6},\n", .{r});
        try writer.print("      \"g_r\": {d:.6},\n", .{g_r});
        try writer.print("      \"h_r\": {d:.6},\n", .{h_r});
        try writer.print("      \"c_r\": {d:.6}\n", .{c_r});

        if (i == oz_solver.grid.n_points - 1) {
            try writer.writeAll("    }\n");
        } else {
            try writer.writeAll("    },\n");
        }
    }
    try writer.writeAll("  ]\n");
    try writer.writeAll("}\n");
}

/// Display correlation functions in formatted console output
pub fn displayRadialFunctions(
    oz_solver: *const Solver,
    max_points: ?usize,
) void {
    const points_to_show = max_points orelse @min(20, oz_solver.grid.n_points);

    std.debug.print("=== Ornstein-Zernike Solution Results ===\n", .{});
    std.debug.print("System: ρ={d:.3}, T={d:.3}, β={d:.3}\n", .{ oz_solver.density, oz_solver.temperature, oz_solver.beta });
    std.debug.print("Grid: {} points, r_max={d:.2}, dr={d:.4}\n", .{ oz_solver.grid.n_points, oz_solver.grid.r_max, oz_solver.grid.dr });
    std.debug.print("\n", .{});

    // Header
    std.debug.print("    r         g_r         h_r         c_r\n", .{});
    std.debug.print("---------------------------------------------\n", .{});

    // Data rows
    for (0..points_to_show) |i| {
        const r = oz_solver.g_r.getRadius(i);
        const g_r = oz_solver.g_r.values[i];
        const h_r = oz_solver.h_r.values[i];
        const c_r = oz_solver.c_r.values[i];

        std.debug.print("{d:8.4}  {d:10.6}  {d:10.6}  {d:10.6}\n", .{ r, g_r, h_r, c_r });
    }

    if (points_to_show < oz_solver.grid.n_points) {
        std.debug.print("... ({} more points)\n", .{oz_solver.grid.n_points - points_to_show});
    }
    std.debug.print("\n", .{});
}

/// Create a simple ASCII plot of g(r)
pub fn plotRadialDistribution(
    oz_solver: *const Solver,
    width: usize,
    height: usize,
) !void {
    const allocator = std.heap.page_allocator;

    // Find min/max values for scaling
    var max_g: f64 = 0.0;
    var min_g: f64 = std.math.inf(f64);

    for (oz_solver.g_r.values) |g_val| {
        if (std.math.isFinite(g_val)) {
            max_g = @max(max_g, g_val);
            min_g = @min(min_g, g_val);
        }
    }

    // Clamp to reasonable range for plotting
    max_g = @min(max_g, 5.0);
    min_g = @max(min_g, 0.0);

    std.debug.print("=== g(r) ASCII Plot ===\n", .{});
    std.debug.print("Range: g(r) ∈ [{d:.2}, {d:.2}]\n", .{ min_g, max_g });
    std.debug.print("\n", .{});

    // Create plot matrix
    const plot = try allocator.alloc([]u8, height);
    defer {
        for (plot) |row| allocator.free(row);
        allocator.free(plot);
    }

    for (0..height) |i| {
        plot[i] = try allocator.alloc(u8, width + 1);
        @memset(plot[i], ' ');
        plot[i][width] = 0; // Null terminator
    }

    // Plot data points
    const step = oz_solver.grid.n_points / width;
    for (0..width) |x| {
        const data_idx = @min(x * step, oz_solver.grid.n_points - 1);
        const g_val = oz_solver.g_r.values[data_idx];

        if (std.math.isFinite(g_val) and g_val >= min_g and g_val <= max_g) {
            const normalized = (g_val - min_g) / (max_g - min_g);
            const y = height - 1 - @as(usize, @intFromFloat(normalized * @as(f64, @floatFromInt(height - 1))));
            plot[y][x] = '*';
        }
    }

    // Print plot with y-axis labels
    for (0..height) |i| {
        const y_val = max_g - (@as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(height - 1))) * (max_g - min_g);
        std.debug.print("{d:5.2} |{s}\n", .{ y_val, plot[i] });
    }

    // Print x-axis
    const dash_line = try allocator.alloc(u8, width + 1);
    defer allocator.free(dash_line);
    @memset(dash_line, '-');
    dash_line[width] = 0;
    std.debug.print("      +{s}\n", .{dash_line});
    std.debug.print("      0{s}{d:.2}\n", .{ "      ", oz_solver.grid.r_max });
    std.debug.print("              r\n", .{});
    std.debug.print("\n", .{});
}

/// Export summary statistics
pub fn exportSummaryStats(
    oz_solver: *const Solver,
    file_path: []const u8,
) !void {
    const file = try std.fs.cwd().createFile(file_path, .{});
    defer file.close();

    const writer = file.writer();

    // Calculate basic statistics
    var g_max: f64 = 0.0;
    var g_min: f64 = std.math.inf(f64);
    var g_mean: f64 = 0.0;
    var finite_count: usize = 0;

    for (oz_solver.g_r.values) |g_val| {
        if (std.math.isFinite(g_val)) {
            g_max = @max(g_max, g_val);
            g_min = @min(g_min, g_val);
            g_mean += g_val;
            finite_count += 1;
        }
    }

    if (finite_count > 0) {
        g_mean /= @as(f64, @floatFromInt(finite_count));
    }

    // Write summary
    try writer.writeAll("Ornstein-Zernike Solution Summary\n");
    try writer.writeAll("=================================\n\n");
    try writer.writeAll("System Parameters:\n");
    try writer.print("  Density (ρ):     {d:.6}\n", .{oz_solver.density});
    try writer.print("  Temperature (T): {d:.6}\n", .{oz_solver.temperature});
    try writer.print("  Beta (β):        {d:.6}\n", .{oz_solver.beta});
    try writer.writeAll("\nGrid Parameters:\n");
    try writer.print("  Points:          {}\n", .{oz_solver.grid.n_points});
    try writer.print("  Maximum r:       {d:.6}\n", .{oz_solver.grid.r_max});
    try writer.print("  Grid spacing:    {d:.6}\n", .{oz_solver.grid.dr});
    try writer.writeAll("\nRadial Distribution Function g(r) Statistics:\n");
    try writer.print("  Minimum:         {d:.6}\n", .{g_min});
    try writer.print("  Maximum:         {d:.6}\n", .{g_max});
    try writer.print("  Mean:            {d:.6}\n", .{g_mean});
    try writer.print("  Finite values:   {}/{}\n", .{ finite_count, oz_solver.grid.n_points });
}
