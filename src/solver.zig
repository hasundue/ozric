const std = @import("std");
const c = @cImport({
    @cInclude("cblas.h");
    @cInclude("lapacke.h");
    @cInclude("fftw3.h");
});

const potentials = @import("potentials.zig");
const grid = @import("grid.zig");
const convergence = @import("convergence.zig");
const closures = @import("closures.zig");
const exact = @import("exact.zig");

pub const Potential = potentials.Potential;
pub const HardSpherePotential = potentials.HardSpherePotential;
pub const LennardJonesPotential = potentials.LennardJonesPotential;
pub const GridParams = grid.GridParams;
pub const RadialFunction = grid.RadialFunction;
pub const ConvergenceParams = convergence.ConvergenceParams;
pub const ConvergenceAccelerator = convergence.ConvergenceAccelerator;

/// Closure relation types for OZ equation
pub const ClosureType = closures.ClosureType;

/// Solver method types
pub const SolverMethod = enum {
    simple_convolution, // Educational simplified convolution
    fft_based, // Proper FFT-based solver
};

pub const Solver = struct {
    allocator: std.mem.Allocator,
    grid: GridParams,

    // Core correlation functions
    h_r: RadialFunction, // h(r) = g(r) - 1, total correlation function
    c_r: RadialFunction, // c(r), direct correlation function
    g_r: RadialFunction, // g(r), radial distribution function

    // FFT workspace for efficient convolution
    fft_workspace: ?[]f64,
    fft_plan_forward: ?c.fftw_plan,
    fft_plan_backward: ?c.fftw_plan,

    // System parameters
    density: f64,
    temperature: f64,
    beta: f64, // β = 1/(kB*T), inverse thermal energy
    potential: ?Potential, // Interaction potential u(r)
    hs_potential_data: ?HardSpherePotential, // Storage for hard sphere potential data
    lj_potential_data: ?LennardJonesPotential, // Storage for LJ potential data

    // Convergence acceleration
    accelerator: ConvergenceAccelerator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, grid_params: GridParams) !Self {
        const h_r = try RadialFunction.init(allocator, grid_params);
        const c_r = try RadialFunction.init(allocator, grid_params);
        const g_r = try RadialFunction.init(allocator, grid_params);

        return Self{
            .allocator = allocator,
            .grid = grid_params,
            .h_r = h_r,
            .c_r = c_r,
            .g_r = g_r,
            .fft_workspace = null,
            .fft_plan_forward = null,
            .fft_plan_backward = null,
            .density = 0.0,
            .temperature = 0.0,
            .beta = 0.0,
            .potential = null,
            .hs_potential_data = null,
            .lj_potential_data = null,
            .accelerator = ConvergenceAccelerator.init(allocator, ConvergenceParams{}),
        };
    }

    pub fn deinit(self: *Self) void {
        self.h_r.deinit();
        self.c_r.deinit();
        self.g_r.deinit();
        self.accelerator.deinit();

        if (self.fft_workspace) |workspace| {
            self.allocator.free(workspace);
        }
        if (self.fft_plan_forward) |plan| {
            c.fftw_destroy_plan(plan);
        }
        if (self.fft_plan_backward) |plan| {
            c.fftw_destroy_plan(plan);
        }
    }

    /// Initialize system with given parameters and potential
    pub fn initSystem(self: *Self, density: f64, temperature: f64, potential: Potential) void {
        self.density = density;
        self.temperature = temperature;
        self.beta = 1.0 / temperature; // Assuming kB = 1 for simplicity
        self.potential = potential;

        // Initialize with a reasonable guess
        self.initMeanFieldGuess();
    }

    /// Initialize system with hard sphere potential
    pub fn initHardSphere(self: *Self, density: f64, temperature: f64, sigma: f64) void {
        // Store the potential data in the solver to ensure lifetime
        self.hs_potential_data = HardSpherePotential.init(sigma);
        const potential = self.hs_potential_data.?.toPotential();
        self.initSystem(density, temperature, potential);

        // Use exact analytical solution as initial guess
        exact.initPYHardSphereGuess(&self.g_r, &self.h_r, &self.c_r, density, sigma);
    }

    /// Initialize system with Lennard-Jones potential
    pub fn initLennardJones(self: *Self, density: f64, temperature: f64, epsilon: f64, sigma: f64) void {
        // Store the potential data in the solver to ensure lifetime
        self.lj_potential_data = LennardJonesPotential.init(epsilon, sigma);
        const potential = self.lj_potential_data.?.toPotential();
        self.initSystem(density, temperature, potential);
    }

    /// Improved initial guess using mean field theory approximation
    pub fn initMeanFieldGuess(self: *Self) void {
        const potential = self.potential orelse return;

        for (0..self.grid.n_points) |i| {
            const r = self.g_r.getRadius(i);
            const u_r = potential.evaluate(r);

            if (std.math.isInf(u_r) and u_r > 0) {
                // Hard core repulsion
                self.g_r.values[i] = 0.0;
            } else if (std.math.isFinite(u_r)) {
                // Mean field approximation: g(r) ≈ exp(-βu_eff(r))
                // where u_eff includes mean field correction
                const mean_field_correction = self.density * 0.1; // Rough approximation
                const u_eff = u_r + mean_field_correction;

                self.g_r.values[i] = @exp(-self.beta * u_eff);

                // Clamp to reasonable values
                if (self.g_r.values[i] > 5.0) self.g_r.values[i] = 5.0;
                if (self.g_r.values[i] < 0.01) self.g_r.values[i] = 0.01;
            } else {
                self.g_r.values[i] = 0.01;
            }

            self.h_r.values[i] = self.g_r.values[i] - 1.0;
        }
    }

    /// Configure convergence acceleration parameters
    pub fn setConvergenceParams(self: *Self, params: ConvergenceParams) void {
        self.accelerator.deinit();
        self.accelerator = ConvergenceAccelerator.init(self.allocator, params);
    }

    /// Initialize FFT plans for FFT-based solver
    pub fn initFFT(self: *Self) !void {
        if (self.fft_workspace != null) return; // Already initialized

        const n = self.grid.n_points;

        // Allocate workspace
        const workspace = try self.allocator.alloc(f64, 2 * n);
        self.fft_workspace = workspace;

        // Create FFT plans for real-to-real transforms
        self.fft_plan_forward = c.fftw_plan_r2r_1d(
            @intCast(n),
            workspace.ptr,
            workspace.ptr,
            c.FFTW_R2HC,
            c.FFTW_ESTIMATE,
        );

        self.fft_plan_backward = c.fftw_plan_r2r_1d(
            @intCast(n),
            workspace.ptr,
            workspace.ptr,
            c.FFTW_HC2R,
            c.FFTW_ESTIMATE,
        );
    }

    /// Main solve method
    pub fn solve(self: *Self, closure: ClosureType, method: SolverMethod, max_iterations: usize, tolerance: f64) !void {
        std.log.info("Starting OZ equation solution with {} closure using {} method...", .{ closure, method });

        // Initialize acceleration methods
        if (method == .fft_based) {
            try self.initFFT();
        }
        try self.accelerator.initAnderson(self.grid.n_points);

        var iteration: usize = 0;
        var err: f64 = std.math.inf(f64);
        var prev_err: f64 = std.math.inf(f64);

        while (iteration < max_iterations and err > tolerance) {
            // Store old h(r) for convergence check
            const old_h = try self.allocator.alloc(f64, self.grid.n_points);
            defer self.allocator.free(old_h);
            @memcpy(old_h, self.h_r.values);

            // Apply closure relation to get c(r)
            try self.applyClosure(closure);

            // Solve OZ equation to get new h(r)
            switch (method) {
                .simple_convolution => try self.solveOZEquationSimple(),
                .fft_based => try self.solveOZEquationFFT(),
            }

            // Apply convergence acceleration
            try self.accelerator.apply(self.h_r.values, old_h, iteration);

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

            // Update adaptive mixing
            self.accelerator.updateAdaptiveMixing(err, prev_err);

            prev_err = err;
            iteration += 1;

            if (iteration % 10 == 0) {
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

        for (0..self.grid.n_points) |i| {
            const r = self.c_r.getRadius(i);
            const u_r = potential.evaluate(r);
            const beta_u = self.beta * u_r;

            // Use the dedicated closure functions
            self.c_r.values[i] = closures.applyClosure(closure, self.g_r.values[i], self.h_r.values[i], beta_u);
        }
    }

    /// Solve OZ equation using simple convolution approximation (educational version)
    fn solveOZEquationSimple(self: *Self) !void {
        // h(r) = c(r) + ρ ∫ c(|r⃗-r⃗'|) h(r') d³r'
        // In spherical coordinates: h(r) = c(r) + 4πρ ∫₀^∞ r'² h(r') [∫₋₁¹ c(s) d(cos θ)] dr'
        // where s = √(r² + r'² - 2rr'cos θ)

        const old_h = try self.allocator.alloc(f64, self.grid.n_points);
        defer self.allocator.free(old_h);
        @memcpy(old_h, self.h_r.values);

        for (0..self.grid.n_points) |i| {
            const r = self.h_r.getRadius(i);
            var convolution: f64 = 0.0;

            // Radial integration over r' with proper 3D geometry
            for (0..self.grid.n_points) |j| {
                const rp = self.h_r.getRadius(j);
                if (rp == 0.0) continue; // Skip r' = 0 point

                // Angular integration: ∫₋₁¹ c(s) d(cos θ) where s = √(r² + r'² - 2rr'cos θ)
                var angular_integral: f64 = 0.0;
                const n_theta = 20; // Number of angular integration points

                for (0..n_theta) |k| {
                    const cos_theta = -1.0 + 2.0 * @as(f64, @floatFromInt(k)) / @as(f64, @floatFromInt(n_theta - 1));

                    // Distance s = |r⃗ - r⃗'|
                    const s_squared = r * r + rp * rp - 2.0 * r * rp * cos_theta;
                    const s = @sqrt(@max(s_squared, 0.0));

                    // Find c(s) by interpolation
                    const s_idx = @min(@as(usize, @intFromFloat(s / self.grid.dr)), self.grid.n_points - 1);
                    const c_s = self.c_r.values[s_idx];

                    angular_integral += c_s * (2.0 / @as(f64, @floatFromInt(n_theta)));
                }

                // Add radial contribution: 4πρ * r'² * h(r') * [angular integral] * dr'
                convolution += 4.0 * std.math.pi * rp * rp * old_h[j] * angular_integral * self.grid.dr;
            }

            self.h_r.values[i] = self.c_r.values[i] + self.density * convolution;
        }
    }

    /// Solve OZ equation using FFT-based convolution (efficient version)
    fn solveOZEquationFFT(self: *Self) !void {
        const workspace = self.fft_workspace orelse return error.NoFFTWorkspace;
        const plan_forward = self.fft_plan_forward orelse return error.NoFFTPlan;
        const plan_backward = self.fft_plan_backward orelse return error.NoFFTPlan;

        const n = self.grid.n_points;

        // Copy c(r) to workspace
        @memcpy(workspace[0..n], self.c_r.values);
        @memset(workspace[n..], 0.0);

        // Transform c(r) to k-space
        c.fftw_execute(plan_forward);

        // Store c(k) temporarily
        const c_k = try self.allocator.alloc(f64, n);
        defer self.allocator.free(c_k);
        @memcpy(c_k, workspace[0..n]);

        // Copy h(r) to workspace and transform
        @memcpy(workspace[0..n], self.h_r.values);
        @memset(workspace[n..], 0.0);
        c.fftw_execute(plan_forward);

        // Multiply in k-space: h(k) = c(k) + ρ * c(k) * h(k)
        // Rearranging: h(k) = c(k) / (1 - ρ * c(k))
        for (0..n) |i| {
            const c_k_val = c_k[i];
            const denominator = 1.0 - self.density * c_k_val;

            if (@abs(denominator) > 1e-10) {
                workspace[i] = c_k_val / denominator;
            } else {
                workspace[i] = c_k_val; // Fallback to prevent division by zero
            }
        }

        // Transform back to real space
        c.fftw_execute(plan_backward);

        // Normalize and copy back
        const norm = 1.0 / @as(f64, @floatFromInt(n));
        for (0..n) |i| {
            self.h_r.values[i] = workspace[i] * norm;
        }
    }
};
