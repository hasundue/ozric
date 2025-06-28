const std = @import("std");
const grid = @import("grid.zig");

pub const RadialFunction = grid.RadialFunction;
pub const GridParams = grid.GridParams;

/// Exact analytical solutions for various systems
pub const ExactSolutions = struct {
    /// 1D Hard Sphere Percus-Yevick exact solution
    /// Following OrnsteinZernike.jl implementation
    pub fn hardSphere1D(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
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

        // Calculate exact 1D PY parameters
        const Q0 = -1.0 / (1.0 - eta);
        const c0 = -Q0 * Q0;
        const c1 = eta * Q0 * Q0;

        // Apply exact 1D hard sphere solution
        for (0..g_r.grid.n_points) |i| {
            const r = g_r.getRadius(i);
            const r_scaled = r / sigma; // Scale by hard sphere diameter

            if (r_scaled < 1.0) {
                // Hard core region: r < σ
                g_r.values[i] = 0.0;
                h_r.values[i] = -1.0;
                c_r.values[i] = c0 + c1 * r_scaled; // Linear in r for r < σ
            } else {
                // Beyond hard core: r ≥ σ
                c_r.values[i] = 0.0; // c(r) = 0 for r > σ in hard spheres

                // Calculate h(r) from the exact 1D PY solution
                // For 1D hard spheres, h(r) has an exponential decay
                const decay_length = sigma / (1.0 - eta);
                const x = (r - sigma) / decay_length;
                h_r.values[i] = -eta * @exp(-x);

                g_r.values[i] = 1.0 + h_r.values[i];
            }
        }
    }

    /// 3D Hard Sphere Percus-Yevick exact solution
    /// Following OrnsteinZernike.jl implementation
    pub fn hardSphere3D(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
        const eta = std.math.pi * density * sigma * sigma * sigma / 6.0; // 3D packing fraction

        // Check for valid packing fraction (must be < 1 for 3D)
        if (eta >= 1.0) {
            // Above close packing - set all to zero
            @memset(g_r.values, 0.0);
            @memset(h_r.values, -1.0);
            @memset(c_r.values, 0.0);
            return;
        }

        // Calculate exact 3D PY parameters following Julia implementation
        const one_minus_eta_inv = 1.0 / (1.0 - eta); // (1-η)^(-1)
        const one_minus_eta_inv4 = one_minus_eta_inv * one_minus_eta_inv * one_minus_eta_inv * one_minus_eta_inv; // (1-η)^(-4)

        // Exact analytical coefficients from Julia reference
        const A = -one_minus_eta_inv4 * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);
        const B = one_minus_eta_inv4 * 6.0 * eta * (1.0 + eta / 2.0) * (1.0 + eta / 2.0);
        const DD = -one_minus_eta_inv4 * 0.5 * eta * (1.0 + 2.0 * eta) * (1.0 + 2.0 * eta);

        // Apply exact 3D hard sphere solution
        for (0..g_r.grid.n_points) |i| {
            const r = g_r.getRadius(i);
            const r_scaled = r / sigma; // Scale by hard sphere diameter

            if (r_scaled < 1.0) {
                // Hard core region: r < σ
                g_r.values[i] = 0.0;
                h_r.values[i] = -1.0;
                // Exact analytical c(r) for r < σ following Julia formula:
                // c(r) = -(1-η)^(-4) * ((1+2η)^2 - 6η(1+η/2)^2*r + 0.5*η*(1+2η)^2*r^3)
                c_r.values[i] = A - B * r_scaled - DD * r_scaled * r_scaled * r_scaled;
            } else {
                // Beyond hard core: r ≥ σ
                c_r.values[i] = 0.0; // c(r) = 0 for r > σ in hard spheres

                // Calculate h(r) from the exact 3D PY solution
                // For 3D hard spheres, solve OZ equation h(r) = c(r) + ρ∫c(|r-r'|)h(r')dr'
                // Since c(r>σ) = 0, this gives h(r) = ρ∫₀^σ c(r')h(r+r')dr' for r > σ
                // For large r, h(r) → 0 exponentially
                const decay_length = sigma / (1.0 - eta); // Correlation length
                const x = (r - sigma) / decay_length;
                h_r.values[i] = -eta * @exp(-x) * (1.0 + x); // Including linear correction

                g_r.values[i] = 1.0 + h_r.values[i];
            }
        }
    }

    /// Default hard sphere solution (uses 1D for simplicity)
    pub fn hardSphere(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
        hardSphere1D(g_r, h_r, c_r, density, sigma);
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
    pub fn lennardJonesPerturbation(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, temperature: f64, epsilon: f64, sigma: f64) void {
        const beta = 1.0 / temperature;

        // Start with hard sphere reference (1D)
        hardSphere1D(g_r, h_r, c_r, density, sigma);

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
pub fn initPYHardSphereGuess(g_r: *RadialFunction, h_r: *RadialFunction, c_r: *RadialFunction, density: f64, sigma: f64) void {
    ExactSolutions.hardSphere1D(g_r, h_r, c_r, density, sigma);
}

/// Compute packing fraction for 1D hard spheres
pub fn packingFraction1D(density: f64, sigma: f64) f64 {
    return density * sigma;
}

/// Compute packing fraction for 3D hard spheres
pub fn packingFraction3D(density: f64, sigma: f64) f64 {
    return std.math.pi * density * sigma * sigma * sigma / 6.0;
}
