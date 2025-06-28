const std = @import("std");

// Re-export all modules for clean public API
pub const potentials = @import("potentials.zig");
pub const grid = @import("grid.zig");
pub const convergence = @import("convergence.zig");
pub const closures = @import("closures.zig");
pub const solver = @import("solver.zig");
pub const export_data = @import("export.zig");

// Re-export commonly used types for convenience
pub const Potential = potentials.Potential;
pub const HardSpherePotential = potentials.HardSpherePotential;
pub const LennardJonesPotential = potentials.LennardJonesPotential;

pub const GridParams = grid.GridParams;
pub const RadialFunction = grid.RadialFunction;
pub const Vector = grid.Vector;

pub const ConvergenceMethod = convergence.ConvergenceMethod;
pub const ConvergenceParams = convergence.ConvergenceParams;
pub const ConvergenceAccelerator = convergence.ConvergenceAccelerator;

pub const ClosureType = solver.ClosureType;
pub const SolverMethod = solver.SolverMethod;
pub const Solver = solver.Solver;

// Keep existing tests for compatibility
test "hard sphere potential with PY closure" {
    const allocator = std.testing.allocator;

    const grid_params = GridParams.init(512, 10.0);
    var oz_solver = try Solver.init(allocator, grid_params);
    defer oz_solver.deinit();

    oz_solver.initHardSphere(0.8, 1.0, 1.0);
    try oz_solver.solve(.percus_yevick, .simple_convolution, 100, 1e-6);
}

test "FFT solver comparison" {
    const allocator = std.testing.allocator;

    const grid_params = GridParams.init(256, 8.0); // Smaller grid for faster test

    // Test simple convolution method
    var solver_simple = try Solver.init(allocator, grid_params);
    defer solver_simple.deinit();

    solver_simple.initHardSphere(0.5, 1.0, 1.0); // Lower density for better convergence
    try solver_simple.solve(.percus_yevick, .simple_convolution, 50, 1e-3);

    // Test FFT-based method
    var solver_fft = try Solver.init(allocator, grid_params);
    defer solver_fft.deinit();

    solver_fft.initHardSphere(0.5, 1.0, 1.0);
    try solver_fft.solve(.percus_yevick, .fft_based, 50, 1e-3);

    // Both methods should produce reasonable results (no crashes, finite values)
    for (0..10) |i| {
        try std.testing.expect(std.math.isFinite(solver_simple.h_r.values[i]));
        try std.testing.expect(std.math.isFinite(solver_fft.h_r.values[i]));
    }
}

test "Lennard-Jones potential with PY closure" {
    const allocator = std.testing.allocator;

    const grid_params = GridParams.init(64, 4.0); // Smaller grid for faster test
    var oz_solver = try Solver.init(allocator, grid_params);
    defer oz_solver.deinit();

    // LJ parameters: moderate conditions for numerical stability
    // density=0.3, temperature=2.0, epsilon=1.0, sigma=1.0
    oz_solver.initLennardJones(0.3, 2.0, 1.0, 1.0);

    // Use conservative convergence parameters
    const stable_params = ConvergenceParams{
        .method = .adaptive_mixing,
        .initial_mixing = 0.05, // Very conservative mixing
        .min_mixing = 0.001,
        .max_mixing = 0.2,
        .error_increase_factor = 1.2,
    };
    oz_solver.setConvergenceParams(stable_params);

    try oz_solver.solve(.percus_yevick, .simple_convolution, 50, 1e-2); // More iterations, looser tolerance

    // Should produce finite results (at least the first few points)
    for (0..3) |i| {
        try std.testing.expect(std.math.isFinite(oz_solver.g_r.values[i]));
    }
}

test "potential evaluation" {
    // Test hard sphere potential
    const hs_pot = HardSpherePotential.init(1.0);
    const hs_potential = hs_pot.toPotential();

    try std.testing.expect(std.math.isInf(hs_potential.evaluate(0.5))); // Inside core
    try std.testing.expectEqual(@as(f64, 0.0), hs_potential.evaluate(1.5)); // Outside core

    // Test Lennard-Jones potential
    const lj_pot = LennardJonesPotential.init(1.0, 1.0);
    const lj_potential = lj_pot.toPotential();

    const u_at_sigma = lj_potential.evaluate(1.0); // At σ
    try std.testing.expectEqual(@as(f64, 0.0), u_at_sigma);

    const u_at_minimum = lj_potential.evaluate(1.122); // At minimum ≈ 2^(1/6)σ
    try std.testing.expect(u_at_minimum < 0.0); // Should be attractive
}

test "convergence acceleration comparison" {
    const allocator = std.testing.allocator;

    const grid_params = GridParams.init(128, 6.0);

    // Test with simple mixing
    var solver_simple = try Solver.init(allocator, grid_params);
    defer solver_simple.deinit();

    const simple_params = ConvergenceParams{
        .method = .simple_mixing,
        .initial_mixing = 0.1,
    };
    solver_simple.setConvergenceParams(simple_params);
    solver_simple.initHardSphere(0.6, 1.0, 1.0); // Slightly lower density for better convergence
    try solver_simple.solve(.percus_yevick, .simple_convolution, 150, 1e-4);

    // Test with adaptive mixing
    var solver_adaptive = try Solver.init(allocator, grid_params);
    defer solver_adaptive.deinit();

    const adaptive_params = ConvergenceParams{
        .method = .adaptive_mixing,
        .initial_mixing = 0.15,
        .min_mixing = 0.01,
        .max_mixing = 0.4,
        .error_increase_factor = 1.5, // More sensitive to error increases
    };
    solver_adaptive.setConvergenceParams(adaptive_params);
    solver_adaptive.initHardSphere(0.6, 1.0, 1.0);
    try solver_adaptive.solve(.percus_yevick, .simple_convolution, 150, 1e-4);

    // Test with Anderson acceleration
    var solver_anderson = try Solver.init(allocator, grid_params);
    defer solver_anderson.deinit();

    const anderson_params = ConvergenceParams{
        .method = .anderson,
        .anderson_depth = 4, // Slightly more history
        .initial_mixing = 0.15,
    };
    solver_anderson.setConvergenceParams(anderson_params);
    solver_anderson.initHardSphere(0.6, 1.0, 1.0);
    try solver_anderson.solve(.percus_yevick, .simple_convolution, 150, 1e-4);

    // All methods should produce finite results
    for (0..5) |i| {
        try std.testing.expect(std.math.isFinite(solver_simple.h_r.values[i]));
        try std.testing.expect(std.math.isFinite(solver_adaptive.h_r.values[i]));
        try std.testing.expect(std.math.isFinite(solver_anderson.h_r.values[i]));
    }
}
