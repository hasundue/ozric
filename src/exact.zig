const std = @import("std");
const grid = @import("grid.zig");
const fourier = @import("fourier.zig");
const oz_fft = @import("oz_fft.zig");

pub const RadialFunction = grid.RadialFunction;
pub const GridParams = grid.GridParams;
const FFTContext = fourier.FFTContext;

/// Exact analytical solutions for various systems using k-space transforms
pub const ExactSolutions = struct {
    /// 1D Hard Sphere Percus-Yevick exact solution using analytical c(k) → g(r) transform
    /// Following OrnsteinZernike.jl implementation with FFT approach
    pub fn hardSphere1D(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) !void {
        const allocator = std.heap.page_allocator; // Use page allocator for simplicity

        // Calculate packing fraction η for 1D hard spheres
        const eta = density * sigma; // 1D packing fraction

        // Check for valid packing fraction (must be < 1 for 1D)
        if (eta >= 1.0) {
            // Above close packing - set all to zero
            @memset(g_r.values, 0.0);
            @memset(h_r.values, -1.0);
            @memset(c_r.values, 0.0);
            return;
        }

        // Step 1: Calculate analytical c(k) in Fourier space
        const n = c_r.grid.n_points;
        const k_max = std.math.pi / c_r.grid.dr; // Nyquist frequency
        const dk = 2.0 * k_max / @as(f64, @floatFromInt(n));

        const c_k = try allocator.alloc(f64, n);
        defer allocator.free(c_k);

        // Exact 1D PY parameters
        const Q0 = -1.0 / (1.0 - eta);
        const c0 = -Q0 * Q0;
        const c1 = eta * Q0 * Q0;

        // Calculate analytical c(k) using exact Fourier transform
        for (0..n) |i| {
            const k = @as(f64, @floatFromInt(i)) * dk;

            if (k < 1e-10) {
                // Handle k → 0 limit
                c_k[i] = c0 + c1 * sigma / 2.0; // Analytical limit
            } else {
                // Exact analytical Fourier transform of 1D PY solution
                const k_sigma = k * sigma;
                const cos_k = @cos(k_sigma);
                const sin_k = @sin(k_sigma);

                c_k[i] = (2.0 * (-c1 + c1 * cos_k + (c0 + c1) * k_sigma * sin_k)) / (k_sigma * k_sigma);
            }
        }

        // Step 2: Use OZ equation to get h(k) from c(k): h(k) = c(k) / (1 - ρ * c(k))
        const h_k = try allocator.alloc(f64, n);
        defer allocator.free(h_k);

        for (0..n) |i| {
            const c_k_val = c_k[i];
            const denominator = 1.0 - density * c_k_val;

            if (@abs(denominator) > 1e-12) {
                h_k[i] = c_k_val / denominator;
            } else {
                h_k[i] = 0.0; // Avoid division by zero
            }
        }

        // Step 3: Inverse FFT to get real-space functions
        var fft_ctx = try FFTContext.init(allocator, n, .estimate);
        defer fft_ctx.deinit();

        // Transform c(k) → c(r)
        try fft_ctx.inverseTransform(c_k, c_r.values);

        // Transform h(k) → h(r)
        try fft_ctx.inverseTransform(h_k, h_r.values);

        // Calculate g(r) = h(r) + 1
        for (0..n) |i| {
            g_r.values[i] = h_r.values[i] + 1.0;
        }

        // Apply hard core constraint: g(r) = 0 for r < σ
        for (0..n) |i| {
            const r = g_r.getRadius(i);
            if (r < sigma) {
                g_r.values[i] = 0.0;
                h_r.values[i] = -1.0;
            }
        }
    }

    /// 3D Hard Sphere Percus-Yevick exact solution using analytical c(k) → g(r) transform
    /// Following OrnsteinZernike.jl implementation with FFT approach
    pub fn hardSphere3D(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) !void {
        const allocator = std.heap.page_allocator; // Use page allocator for simplicity

        const eta = std.math.pi * density * sigma * sigma * sigma / 6.0; // 3D packing fraction

        // Check for valid packing fraction (must be < 1 for 3D)
        if (eta >= 1.0) {
            // Above close packing - set all to zero
            @memset(g_r.values, 0.0);
            @memset(h_r.values, -1.0);
            @memset(c_r.values, 0.0);
            return;
        }

        // Step 1: Calculate analytical c(k) in Fourier space
        const n = c_r.grid.n_points;
        const k_max = std.math.pi / c_r.grid.dr; // Nyquist frequency
        const dk = 2.0 * k_max / @as(f64, @floatFromInt(n));

        const c_k = try allocator.alloc(f64, n);
        defer allocator.free(c_k);

        // Calculate exact 3D PY parameters following Julia implementation
        const one_minus_eta_inv4 = std.math.pow(f64, 1.0 - eta, -4.0); // (1-η)^(-4)

        // Exact analytical coefficients from Julia reference
        const A = -one_minus_eta_inv4 * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);
        const B = one_minus_eta_inv4 * 6.0 * eta * (1.0 + eta / 2.0) * (1.0 + eta / 2.0);
        const DD = -one_minus_eta_inv4 * 0.5 * eta * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);

        // Calculate analytical c(k) using exact 3D Fourier transform from Julia
        for (0..n) |i| {
            const k = @as(f64, @floatFromInt(i)) * dk;

            if (k < 1e-10) {
                // Handle k → 0 limit using L'Hôpital's rule
                c_k[i] = (4.0 * std.math.pi / 6.0) * (24.0 * DD + 2.0 * B + A); // Analytical limit
            } else {
                // Exact analytical 3D Fourier transform
                const k_sigma = k * sigma;
                const k2 = k * k;
                const k4 = k2 * k2;
                const k6 = k4 * k2;
                const cos_k = @cos(k_sigma);
                const sin_k = @sin(k_sigma);

                const term1 = 24.0 * DD - 2.0 * B * k2;
                const term2 = -(24.0 * DD - 2.0 * (B + 6.0 * DD) * k2 + (A + B + DD) * k4) * cos_k;
                const term3 = k_sigma * (-24.0 * DD + (A + 2.0 * B + 4.0 * DD) * k2) * sin_k;

                c_k[i] = 4.0 * std.math.pi * (term1 + term2 + term3) / k6;
            }
        }

        // Step 2: Use OZ equation to get h(k) from c(k): h(k) = c(k) / (1 - ρ * c(k))
        const h_k = try allocator.alloc(f64, n);
        defer allocator.free(h_k);

        for (0..n) |i| {
            const c_k_val = c_k[i];
            const denominator = 1.0 - density * c_k_val;

            if (@abs(denominator) > 1e-12) {
                h_k[i] = c_k_val / denominator;
            } else {
                h_k[i] = 0.0; // Avoid division by zero
            }
        }

        // Step 3: Inverse FFT to get real-space functions
        var fft_ctx = try FFTContext.init(allocator, n, .estimate);
        defer fft_ctx.deinit();

        // Transform c(k) → c(r)
        try fft_ctx.inverseTransform(c_k, c_r.values);

        // Transform h(k) → h(r)
        try fft_ctx.inverseTransform(h_k, h_r.values);

        // Calculate g(r) = h(r) + 1
        for (0..n) |i| {
            g_r.values[i] = h_r.values[i] + 1.0;
        }

        // Apply hard core constraint: g(r) = 0 for r < σ
        for (0..n) |i| {
            const r = g_r.getRadius(i);
            if (r < sigma) {
                g_r.values[i] = 0.0;
                h_r.values[i] = -1.0;
            }
        }
    }

    /// Default hard sphere solution (uses 1D for simplicity)
    pub fn hardSphere(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) !void {
        try hardSphere1D(g_r, h_r, c_r, density, sigma);
    }

    /// Ideal gas exact solution (trivial case)
    pub fn idealGas(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction) void {
        for (0..g_r.grid.n_points) |i| {
            g_r.values[i] = 1.0; // No correlations
            h_r.values[i] = 0.0; // h(r) = g(r) - 1
            c_r.values[i] = 0.0; // No direct correlations
        }
    }

    /// Calculate 1D hard sphere direct correlation function in k-space
    /// Following the exact analytical Fourier transform
    pub fn hardSphere1DFourier(c_k: []f64, k_values: []const f64, density: f64, sigma: f64) void {
        const eta = density * sigma; // 1D packing fraction

        if (eta >= 1.0) {
            @memset(c_k, 0.0);
            return;
        }

        // Exact 1D PY parameters
        const Q0 = -1.0 / (1.0 - eta);
        const c0 = -Q0 * Q0;
        const c1 = eta * Q0 * Q0;

        // Apply exact Fourier transform of c(r)
        for (0..k_values.len) |i| {
            const k = k_values[i];

            if (k < 1e-10) {
                // Handle k → 0 limit
                c_k[i] = c0 + c1 / 2.0; // Analytical limit
            } else {
                // Exact analytical Fourier transform
                const k_sigma = k * sigma;
                const cos_k = @cos(k_sigma);
                const sin_k = @sin(k_sigma);

                c_k[i] = (2.0 * (-c1 + c1 * cos_k + (c0 + c1) * k_sigma * sin_k)) / (k_sigma * k_sigma);
            }
        }
    }

    /// Calculate 3D hard sphere direct correlation function in k-space
    /// Following the exact analytical Fourier transform
    pub fn hardSphere3DFourier(c_k: []f64, k_values: []const f64, density: f64, sigma: f64) void {
        const eta = std.math.pi * density * sigma * sigma * sigma / 6.0; // 3D packing fraction

        if (eta >= 1.0) {
            @memset(c_k, 0.0);
            return;
        }

        // Calculate exact 3D PY parameters following Julia implementation
        const one_minus_eta_inv4 = 1.0 / std.math.pow(f64, 1.0 - eta, 4.0); // (1-η)^(-4)

        // Exact analytical coefficients from Julia reference
        const A = -one_minus_eta_inv4 * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);
        const B = one_minus_eta_inv4 * 6.0 * eta * (1.0 + eta / 2.0) * (1.0 + eta / 2.0);
        const DD = -one_minus_eta_inv4 * 0.5 * eta * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);

        // Apply exact Fourier transform of c(r)
        for (0..k_values.len) |i| {
            const k = k_values[i];

            if (k < 1e-10) {
                // Handle k → 0 limit
                c_k[i] = A + B / 2.0 + DD / 4.0; // Analytical limit
            } else {
                // Exact analytical Fourier transform from Julia reference:
                // Ck = 4π/k^6 * (24*DD - 2*B*k^2 - (24*DD - 2*(B+6*DD)*k^2 + (A+B+DD)*k^4)*cos(k) + k*(-24*DD + (A+2*B+4*DD)*k^2)*sin(k))
                const k_sigma = k * sigma;
                const k2 = k_sigma * k_sigma;
                const k4 = k2 * k2;
                const k6 = k4 * k2;
                const cos_k = @cos(k_sigma);
                const sin_k = @sin(k_sigma);

                const term1 = 24.0 * DD - 2.0 * B * k2;
                const term2 = -(24.0 * DD - 2.0 * (B + 6.0 * DD) * k2 + (A + B + DD) * k4) * cos_k;
                const term3 = k_sigma * (-24.0 * DD + (A + 2.0 * B + 4.0 * DD) * k2) * sin_k;

                c_k[i] = 4.0 * std.math.pi / k6 * (term1 + term2 + term3);
            }
        }
    }

    /// Lennard-Jones approximate solution using perturbation theory
    pub fn lennardJonesPerturbation(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, temperature: f64, epsilon: f64, sigma: f64) !void {
        const beta = 1.0 / temperature;

        // Start with hard sphere reference (1D)
        try hardSphere1D(g_r, h_r, c_r, density, sigma);

        // Apply first-order perturbation correction
        for (0..g_r.grid.n_points) |i| {
            const r = g_r.getRadius(i);

            if (r >= sigma) {
                // LJ attractive tail: u_attr(r) = -4ε[(σ/r)^6]
                const sigma_over_r = sigma / r;
                const sigma_over_r_6 = std.math.pow(f64, sigma_over_r, 6.0);
                const u_attr = -4.0 * epsilon * sigma_over_r_6;

                // First-order perturbation: g(r) ≈ g_hs(r) * exp(-βu_attr(r))
                const perturbation = @exp(-beta * u_attr);

                g_r.values[i] *= perturbation;
                h_r.values[i] = g_r.values[i] - 1.0;

                // Approximate c(r) correction (simplified)
                if (r < 3.0 * sigma) {
                    c_r.values[i] += -beta * u_attr * 0.1; // Rough approximation
                }
            }
        }
    }
};

/// Initialize with 1D Percus-Yevick exact solution for hard spheres
pub fn initPYHardSphereGuess(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) !void {
    try ExactSolutions.hardSphere1D(g_r, h_r, c_r, density, sigma);
}

/// Compute packing fraction for 1D hard spheres
pub fn packingFraction1D(density: f64, sigma: f64) f64 {
    return density * sigma;
}

/// Compute packing fraction for 3D hard spheres
pub fn packingFraction3D(density: f64, sigma: f64) f64 {
    return std.math.pi * density * sigma * sigma * sigma / 6.0;
}
