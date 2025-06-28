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

    // System parameters
    density: f64,
    temperature: f64,

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
            .density = 0.0,
            .temperature = 0.0,
        };
    }

    pub fn deinit(self: *Self) void {
        self.allocator.free(self.workspace);
        self.h_r.deinit();
        self.c_r.deinit();
        self.g_r.deinit();
    }

    /// Initialize with hard sphere potential as starting guess
    pub fn initHardSphere(self: *Self, density: f64, temperature: f64, sigma: f64) void {
        self.density = density;
        self.temperature = temperature;

        // Initialize g(r) with hard sphere step function
        for (0..self.grid.n_points) |i| {
            const r = self.g_r.getRadius(i);
            if (r < sigma) {
                self.g_r.values[i] = 0.0;
            } else {
                self.g_r.values[i] = 1.0;
            }
            self.h_r.values[i] = self.g_r.values[i] - 1.0;
        }
    }

    /// Solve the Ornstein-Zernike equation using specified closure relation
    pub fn solve(self: *Self, closure: ClosureType, max_iterations: usize, tolerance: f64) !void {
        std.log.info("Starting OZ equation solution with {} closure...", .{closure});

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
            try self.solveOZEquation();

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
        switch (closure) {
            .hypernetted_chain => {
                // HNC: c(r) = h(r) - ln(g(r))
                for (0..self.grid.n_points) |i| {
                    if (self.g_r.values[i] > 0.0) {
                        self.c_r.values[i] = self.h_r.values[i] - @log(self.g_r.values[i]);
                    } else {
                        self.c_r.values[i] = -std.math.inf(f64);
                    }
                }
            },
            .percus_yevick => {
                // PY: c(r) = (1 - 1/g(r)) * h(r) for hard spheres
                for (0..self.grid.n_points) |i| {
                    if (self.g_r.values[i] > 0.0) {
                        self.c_r.values[i] = (1.0 - 1.0 / self.g_r.values[i]) * self.h_r.values[i];
                    } else {
                        self.c_r.values[i] = 0.0;
                    }
                }
            },
        }
    }

    /// Solve OZ equation using simple convolution approximation for now
    fn solveOZEquation(self: *Self) !void {
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
    try solver.solve(.percus_yevick, 10, 1e-6);
}
