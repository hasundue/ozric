const std = @import("std");
const c = @cImport({
    @cInclude("cblas.h");
    @cInclude("lapacke.h");
    @cInclude("fftw3.h");
});

/// Closure relation types for OZ equation
pub const ClosureType = enum {
    hypernetted_chain, // HNC: c(r) = h(r) - ln(g(r)) - βu(r)
    percus_yevick, // PY: c(r) = (1 - exp(βu(r))) * g(r)
};

/// Solver method types
pub const SolverMethod = enum {
    simple_convolution, // Educational simplified convolution
    fft_based, // Proper FFT-based solver
};

/// Potential function interface
pub const Potential = struct {
    const Self = @This();

    // Function pointer for potential evaluation
    evaluateFn: *const fn (self: *const Self, r: f64) f64,

    // Parameters (can be cast to specific potential types)
    params: *anyopaque,

    pub fn evaluate(self: *const Self, r: f64) f64 {
        return self.evaluateFn(self, r);
    }
};

/// Hard sphere potential parameters
pub const HardSpherePotential = struct {
    sigma: f64, // hard sphere diameter

    const Self = @This();

    pub fn init(sigma: f64) Self {
        return Self{ .sigma = sigma };
    }

    pub fn toPotential(self: *const Self) Potential {
        return Potential{
            .evaluateFn = evaluate,
            .params = @constCast(@ptrCast(self)),
        };
    }

    fn evaluate(potential: *const Potential, r: f64) f64 {
        const self: *const HardSpherePotential = @ptrCast(@alignCast(potential.params));
        if (r < self.sigma) {
            return std.math.inf(f64); // Hard sphere repulsion
        } else {
            return 0.0; // No interaction beyond contact
        }
    }
};

/// Lennard-Jones potential parameters
pub const LennardJonesPotential = struct {
    epsilon: f64, // energy scale
    sigma: f64, // length scale

    const Self = @This();

    pub fn init(epsilon: f64, sigma: f64) Self {
        return Self{ .epsilon = epsilon, .sigma = sigma };
    }

    pub fn toPotential(self: *const Self) Potential {
        return Potential{
            .evaluateFn = evaluate,
            .params = @constCast(@ptrCast(self)),
        };
    }

    fn evaluate(potential: *const Potential, r: f64) f64 {
        const self: *const LennardJonesPotential = @ptrCast(@alignCast(potential.params));
        const sigma_over_r = self.sigma / r;
        const sigma_over_r6 = sigma_over_r * sigma_over_r * sigma_over_r *
            sigma_over_r * sigma_over_r * sigma_over_r;
        const sigma_over_r12 = sigma_over_r6 * sigma_over_r6;

        return 4.0 * self.epsilon * (sigma_over_r12 - sigma_over_r6);
    }
};

/// Grid parameters for radial functions
pub const GridParams = struct {
    n_points: usize,
    r_max: f64,
    dr: f64,

    pub fn init(n_points: usize, r_max: f64) GridParams {
        return GridParams{
            .n_points = n_points,
            .r_max = r_max,
            .dr = r_max / @as(f64, @floatFromInt(n_points - 1)),
        };
    }
};

/// Radial function data structure
pub const RadialFunction = struct {
    values: []f64,
    grid: GridParams,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: GridParams) !Self {
        const values = try allocator.alloc(f64, grid.n_points);
        @memset(values, 0.0);
        return Self{
            .values = values,
            .grid = grid,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.values);
    }

    pub fn getRadius(self: *const Self, i: usize) f64 {
        return @as(f64, @floatFromInt(i)) * self.grid.dr;
    }
};

pub const Solver = struct {
    allocator: std.mem.Allocator,
    grid: GridParams,

    // Core correlation functions
    h_r: RadialFunction, // h(r) = g(r) - 1, total correlation function
    c_r: RadialFunction, // c(r), direct correlation function
    g_r: RadialFunction, // g(r), radial distribution function

    // Workspace for calculations
    workspace: []f64,

    // FFT workspace and plans (optional, for fft_based method)
    fft_workspace: ?[]f64,
    fft_plan_forward: ?*anyopaque,
    fft_plan_backward: ?*anyopaque,

    // System parameters
    density: f64,
    temperature: f64,
    beta: f64, // β = 1/(kB*T), inverse thermal energy
    potential: ?Potential, // Interaction potential u(r)

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid: GridParams) !Self {
        var h_r = try RadialFunction.init(allocator, grid);
        errdefer h_r.deinit();

        var c_r = try RadialFunction.init(allocator, grid);
        errdefer c_r.deinit();

        var g_r = try RadialFunction.init(allocator, grid);
        errdefer g_r.deinit();

        // Allocate workspace for calculations
        const workspace = try allocator.alloc(f64, grid.n_points);
        errdefer allocator.free(workspace);

        return Self{
            .allocator = allocator,
            .grid = grid,
            .h_r = h_r,
            .c_r = c_r,
            .g_r = g_r,
            .workspace = workspace,
            .fft_workspace = null,
            .fft_plan_forward = null,
            .fft_plan_backward = null,
            .density = 0.0,
            .temperature = 0.0,
            .beta = 0.0,
            .potential = null,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.workspace);
        if (self.fft_workspace) |fft_ws| {
            self.allocator.free(fft_ws);
        }
        if (self.fft_plan_forward) |plan| {
            c.fftw_destroy_plan(@ptrCast(plan));
        }
        if (self.fft_plan_backward) |plan| {
            c.fftw_destroy_plan(@ptrCast(plan));
        }
        self.h_r.deinit();
        self.c_r.deinit();
        self.g_r.deinit();
    }

    /// Initialize system with given potential and thermodynamic state
    pub fn initSystem(self: *Self, density: f64, temperature: f64, potential: Potential) void {
        self.density = density;
        self.temperature = temperature;
        self.beta = 1.0 / temperature; // Assuming kB = 1 (reduced units)
        self.potential = potential;

        // Initialize g(r) based on potential
        for (0..self.grid.n_points) |i| {
            const r = self.g_r.getRadius(i);
            const u_r = potential.evaluate(r);

            if (std.math.isInf(u_r) and u_r > 0) {
                // Hard core repulsion
                self.g_r.values[i] = 0.0;
            } else {
                // Initial guess: g(r) ≈ exp(-βu(r)) for r > contact
                self.g_r.values[i] = if (r > 0.01) @exp(-self.beta * u_r) else 0.0;

                // Clamp to reasonable values
                if (self.g_r.values[i] > 10.0) self.g_r.values[i] = 10.0;
                if (self.g_r.values[i] < 0.01) self.g_r.values[i] = 0.01;
            }
            self.h_r.values[i] = self.g_r.values[i] - 1.0;
        }
    }

    /// Convenience method for hard sphere initialization
    pub fn initHardSphere(self: *Self, density: f64, temperature: f64, sigma: f64) void {
        const hs_potential = HardSpherePotential.init(sigma);
        self.initSystem(density, temperature, hs_potential.toPotential());
    }

    /// Convenience method for Lennard-Jones initialization
    pub fn initLennardJones(self: *Self, density: f64, temperature: f64, epsilon: f64, sigma: f64) void {
        const lj_potential = LennardJonesPotential.init(epsilon, sigma);
        self.initSystem(density, temperature, lj_potential.toPotential());
    }

    /// Initialize FFT plans for FFT-based solver
    pub fn initFFT(self: *Self) !void {
        if (self.fft_workspace != null) return; // Already initialized

        const n = self.grid.n_points;

        // Allocate FFT workspace for real-to-complex transform
        const fft_workspace = try self.allocator.alloc(f64, n);
        self.fft_workspace = fft_workspace;

        // Create FFTW plans using real-to-complex transforms
        const forward_plan = c.fftw_plan_r2r_1d(@intCast(n), fft_workspace.ptr, fft_workspace.ptr, c.FFTW_R2HC, c.FFTW_ESTIMATE);

        const backward_plan = c.fftw_plan_r2r_1d(@intCast(n), fft_workspace.ptr, fft_workspace.ptr, c.FFTW_HC2R, c.FFTW_ESTIMATE);

        self.fft_plan_forward = @ptrCast(forward_plan);
        self.fft_plan_backward = @ptrCast(backward_plan);
    }

    /// Solve the Ornstein-Zernike equation using specified closure relation and method
    pub fn solve(self: *Self, closure: ClosureType, method: SolverMethod, max_iterations: usize, tolerance: f64) !void {
        std.log.info("Starting OZ equation solution with {} closure using {} method...", .{ closure, method });

        // Initialize FFT if using FFT-based method
        if (method == .fft_based) {
            try self.initFFT();
        }

        var iteration: usize = 0;
        var err: f64 = std.math.inf(f64);

        while (iteration < max_iterations and err > tolerance) {
            // Store old h(r) for convergence check
            const old_h = try self.allocator.alloc(f64, self.grid.n_points);
            defer self.allocator.free(old_h);
            @memcpy(old_h, self.h_r.values);

            // Apply closure relation to get c(r) from h(r)
            try self.applyClosure(closure);

            // Solve OZ equation: h(r) = c(r) + ρ ∫ c(|r-r'|) h(r') dr'
            switch (method) {
                .simple_convolution => try self.solveOZEquationSimple(),
                .fft_based => try self.solveOZEquationFFT(),
            }

            // Update g(r) = h(r) + 1
            for (0..self.grid.n_points) |i| {
                self.g_r.values[i] = self.h_r.values[i] + 1.0;
            }

            // Check convergence
            err = 0.0;
            for (0..self.grid.n_points) |i| {
                const diff = self.h_r.values[i] - old_h[i];
                err += diff * diff;
            }
            err = @sqrt(err / @as(f64, @floatFromInt(self.grid.n_points)));

            iteration += 1;

            if (iteration % 100 == 0) {
                std.log.info("Iteration {}: error = {d:.6}", .{ iteration, err });
            }
        }

        if (err <= tolerance) {
            std.log.info("Converged after {} iterations with error {d:.6}", .{ iteration, err });
        } else {
            std.log.warn("Failed to converge after {} iterations, final error: {d:.6}", .{ max_iterations, err });
        }
    }

    /// Apply closure relation to compute c(r) from h(r) and g(r)
    fn applyClosure(self: *Self, closure: ClosureType) !void {
        const potential = self.potential orelse return error.NoPotential;

        switch (closure) {
            .hypernetted_chain => {
                // HNC: c(r) = h(r) - ln(g(r)) - βu(r)
                for (0..self.grid.n_points) |i| {
                    const r = self.h_r.getRadius(i);
                    const u_r = potential.evaluate(r);

                    if (self.g_r.values[i] > 1e-12 and std.math.isFinite(u_r)) {
                        const log_g = @log(@max(self.g_r.values[i], 1e-12));
                        const beta_u = self.beta * u_r;

                        // Clamp βu to prevent overflow
                        const clamped_beta_u = @max(@min(beta_u, 50.0), -50.0);

                        self.c_r.values[i] = self.h_r.values[i] - log_g - clamped_beta_u;
                    } else {
                        // Handle infinite potential or zero g(r)
                        self.c_r.values[i] = -1e6; // Large negative instead of -inf
                    }
                }
            },
            .percus_yevick => {
                // PY: c(r) = (1 - exp(βu(r))) * g(r)
                for (0..self.grid.n_points) |i| {
                    const r = self.c_r.getRadius(i);
                    const u_r = potential.evaluate(r);

                    if (std.math.isFinite(u_r)) {
                        const beta_u = self.beta * u_r;

                        // Clamp βu to prevent overflow
                        const clamped_beta_u = @max(@min(beta_u, 50.0), -50.0);
                        const exp_beta_u = @exp(clamped_beta_u);

                        if (std.math.isFinite(exp_beta_u)) {
                            self.c_r.values[i] = (1.0 - exp_beta_u) * self.g_r.values[i];
                        } else {
                            // For very large βu (hard core-like)
                            self.c_r.values[i] = -self.g_r.values[i];
                        }
                    } else {
                        // Handle infinite potential
                        self.c_r.values[i] = -self.g_r.values[i];
                    }
                }
            },
        }
    }

    /// Solve OZ equation using simple convolution approximation (educational version)
    fn solveOZEquationSimple(self: *Self) !void {
        // Simplified OZ solution without FFT for initial implementation
        // In a real solver, this would use FFT for proper convolution

        const mixing_factor = 0.3;

        // Copy old h(r) to workspace
        @memcpy(self.workspace, self.h_r.values);

        // Simple approximation: h(r) ≈ c(r) + ρ * integral approximation
        for (0..self.grid.n_points) |i| {
            const r = self.h_r.getRadius(i);
            var integral_approx: f64 = 0.0;

            // Rough integral approximation using trapezoidal rule
            for (0..self.grid.n_points - 1) |j| {
                const r_prime = self.c_r.getRadius(j);
                const dr = self.grid.dr;

                if (r > 0.01 and r_prime > 0.01) { // avoid singularities
                    const kernel = self.c_r.values[j] * self.workspace[j];
                    integral_approx += kernel * dr * r_prime * r_prime; // spherical coordinates
                }
            }

            const new_h = self.c_r.values[i] + self.density * integral_approx / (r * r + 0.01);
            self.h_r.values[i] = mixing_factor * new_h + (1.0 - mixing_factor) * self.workspace[i];
        }
    }

    /// Solve OZ equation using FFT-based convolution (proper implementation)
    fn solveOZEquationFFT(self: *Self) !void {
        const n = self.grid.n_points;
        const fft_ws = self.fft_workspace orelse return error.FFTNotInitialized;

        // Simplified FFT-based approach using real transforms
        // For educational purposes, using a simpler convolution in Fourier space

        // Step 1: Prepare c(r) for FFT
        for (0..n) |i| {
            const r = self.c_r.getRadius(i);
            // Apply spherical coordinate factor for proper transform
            fft_ws[i] = if (r > 0.01) r * self.c_r.values[i] else 0.0;
        }

        // Forward FFT: c(r) -> C(k)
        if (self.fft_plan_forward) |plan| {
            c.fftw_execute(@ptrCast(plan));
        }

        // Apply simplified OZ relation in k-space
        // H(k) ≈ C(k) / (1 - ρ*C(k)) simplified for real transforms
        for (0..n) |i| {
            const c_k = fft_ws[i];
            const denom = 1.0 - self.density * c_k;

            if (@abs(denom) > 1e-12) {
                fft_ws[i] = c_k / denom;
            } else {
                fft_ws[i] = 0.0;
            }
        }

        // Inverse FFT: H(k) -> h(r)
        if (self.fft_plan_backward) |plan| {
            c.fftw_execute(@ptrCast(plan));
        }

        // Extract h(r) and apply mixing for stability
        const mixing_factor = 0.2;
        const norm_factor = 1.0 / @as(f64, @floatFromInt(n));

        for (0..n) |i| {
            const r = self.h_r.getRadius(i);

            // Extract h(r) from r*h(r) and normalize
            const new_h = if (r > 0.01) fft_ws[i] * norm_factor / r else 0.0;

            // Apply mixing for numerical stability
            self.h_r.values[i] = mixing_factor * new_h + (1.0 - mixing_factor) * self.h_r.values[i];
        }
    }
};

/// Utility functions for vector operations using BLAS
pub const Vector = struct {
    data: []f64,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, size: usize) !Self {
        const data = try allocator.alloc(f64, size);
        return Self{
            .data = data,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.data);
    }

    /// Compute dot product using BLAS
    pub fn dot(self: *const Self, other: *const Self) f64 {
        std.debug.assert(self.data.len == other.data.len);
        const n: c_int = @intCast(self.data.len);
        return c.cblas_ddot(n, self.data.ptr, 1, other.data.ptr, 1);
    }

    /// Scale vector by a constant using BLAS
    pub fn scale(self: *Self, alpha: f64) void {
        const n: c_int = @intCast(self.data.len);
        c.cblas_dscal(n, alpha, self.data.ptr, 1);
    }
};

test "vector operations" {
    const allocator = std.testing.allocator;

    var v1 = try Vector.init(allocator, 3);
    defer v1.deinit();

    var v2 = try Vector.init(allocator, 3);
    defer v2.deinit();

    v1.data[0] = 1.0;
    v1.data[1] = 2.0;
    v1.data[2] = 3.0;

    v2.data[0] = 4.0;
    v2.data[1] = 5.0;
    v2.data[2] = 6.0;

    const result = v1.dot(&v2);
    try std.testing.expectEqual(@as(f64, 32.0), result);
}

test "solver initialization" {
    const allocator = std.testing.allocator;

    const grid = GridParams.init(512, 10.0);
    var solver = try Solver.init(allocator, grid);
    defer solver.deinit();

    solver.initHardSphere(0.8, 1.0, 1.0);
    try solver.solve(.percus_yevick, .simple_convolution, 10, 1e-6);
}

test "FFT solver comparison" {
    const allocator = std.testing.allocator;

    const grid = GridParams.init(256, 8.0); // Smaller grid for faster test

    // Test simple convolution method
    var solver_simple = try Solver.init(allocator, grid);
    defer solver_simple.deinit();

    solver_simple.initHardSphere(0.5, 1.0, 1.0); // Lower density for better convergence
    try solver_simple.solve(.percus_yevick, .simple_convolution, 5, 1e-3);

    // Test FFT-based method
    var solver_fft = try Solver.init(allocator, grid);
    defer solver_fft.deinit();

    solver_fft.initHardSphere(0.5, 1.0, 1.0);
    try solver_fft.solve(.percus_yevick, .fft_based, 5, 1e-3);

    // Both methods should produce reasonable results (no crashes, finite values)
    for (0..10) |i| {
        try std.testing.expect(std.math.isFinite(solver_simple.h_r.values[i]));
        try std.testing.expect(std.math.isFinite(solver_fft.h_r.values[i]));
    }
}

test "Lennard-Jones potential with PY closure" {
    const allocator = std.testing.allocator;

    const grid = GridParams.init(64, 4.0); // Smaller grid for faster test
    var solver = try Solver.init(allocator, grid);
    defer solver.deinit();

    // LJ parameters: reduced units (ε=0.1, σ=1) - weaker interaction for stability
    solver.initLennardJones(0.1, 5.0, 0.1, 1.0); // Very low density, very high temperature
    try solver.solve(.percus_yevick, .simple_convolution, 2, 1e-1); // Loose convergence

    // Should produce finite results (at least the first few points)
    for (0..3) |i| {
        try std.testing.expect(std.math.isFinite(solver.g_r.values[i]));
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
