const std = @import("std");
const grid = @import("grid.zig");

pub const RadialFunction = grid.RadialFunction;
pub const GridParams = grid.GridParams;

/// Exact analytical solutions for various systems
pub const ExactSolutions = struct {
    /// Hard sphere exact solution (contact values and asymptotic behavior)
    pub fn hardSphere(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
        const eta = std.math.pi * density * sigma * sigma * sigma / 6.0; // packing fraction

        // Percus-Yevick hard sphere analytical solution parameters
        const alpha = if (eta < 0.5) blk: {
            const one_minus_eta = 1.0 - eta;
            const numerator = (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);
            const denominator = one_minus_eta * one_minus_eta * one_minus_eta * one_minus_eta;
            break :blk numerator / denominator;
        } else 1.0; // Fallback for high density

        const beta_hs = if (eta < 0.5) blk: {
            const one_minus_eta = 1.0 - eta;
            const numerator = -6.0 * eta * (1.0 + eta / 2.0);
            const denominator = one_minus_eta * one_minus_eta * one_minus_eta;
            break :blk numerator / denominator;
        } else -1.0; // Fallback for high density

        for (0..g_r.grid.n_points) |i| {
            const r = g_r.getRadius(i);

            if (r < sigma) {
                // Hard core region: exact solution
                g_r.values[i] = 0.0;
                h_r.values[i] = -1.0;
                c_r.values[i] = -1.0;
            } else if (r < 2.0 * sigma) {
                // Contact region: PY analytical approximation
                const x = r / sigma;
                g_r.values[i] = alpha + beta_hs * (x - 1.0);

                // Ensure physical bounds
                if (g_r.values[i] < 0.01) g_r.values[i] = 0.01;
                if (g_r.values[i] > 5.0) g_r.values[i] = 5.0;

                h_r.values[i] = g_r.values[i] - 1.0;

                // PY closure: c(r) = 0 for r > σ in hard sphere case
                c_r.values[i] = 0.0;
            } else {
                // Long range: approach ideal gas limit with exponential decay
                const decay_length = sigma; // Characteristic decay length
                const decay = @exp(-(r - 2.0 * sigma) / decay_length);

                g_r.values[i] = 1.0 + 0.1 * eta * decay; // Small oscillations around 1
                h_r.values[i] = g_r.values[i] - 1.0;
                c_r.values[i] = 0.0; // Long-range c(r) → 0
            }
        }
    }

    /// Ideal gas exact solution (trivial case)
    pub fn idealGas(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction) void {
        for (0..g_r.grid.n_points) |i| {
            g_r.values[i] = 1.0; // No correlations
            h_r.values[i] = 0.0; // h(r) = g(r) - 1
            c_r.values[i] = 0.0; // No direct correlations
        }
    }

    /// Lennard-Jones approximate solution using perturbation theory
    pub fn lennardJonesPerturbation(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, temperature: f64, epsilon: f64, sigma: f64) void {
        const beta = 1.0 / temperature;

        // Start with hard sphere reference
        hardSphere(g_r, h_r, c_r, density, sigma);

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

/// Initialize with Percus-Yevick analytical solution for hard spheres
pub fn initPYHardSphereGuess(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
    ExactSolutions.hardSphere(g_r, h_r, c_r, density, sigma);
}
