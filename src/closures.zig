const std = @import("std");

/// Closure relation types for Ornstein-Zernike equations
pub const ClosureType = enum {
    hypernetted_chain, // HNC
    percus_yevick, // PY
};

/// Apply Hypernetted Chain (HNC) closure relation
/// HNC: c(r) = h(r) - ln(g(r)) - βu(r)
pub fn applyHNC(g_r: f64, h_r: f64, beta_u: f64) f64 {
    if (g_r > 0.0) {
        if (std.math.isInf(beta_u)) {
            // For infinite potentials, HNC gives specific behavior
            if (beta_u > 0) {
                // Repulsive infinite potential: c(r) = h(r) - ln(g(r)) - ∞
                return h_r - @log(g_r) - 1e10;
            } else {
                // Attractive infinite potential (shouldn't happen)
                return h_r - @log(g_r) + 1e10;
            }
        } else {
            // Normal case with finite potential
            return h_r - @log(g_r) - beta_u;
        }
    } else {
        // g(r) ≤ 0: HNC becomes problematic, use limiting case
        return h_r;
    }
}

/// Apply Percus-Yevick (PY) closure relation
/// PY: c(r) = (1 - exp(βu(r))) * g(r)
pub fn applyPY(g_r: f64, beta_u: f64) f64 {
    if (std.math.isInf(beta_u)) {
        // For infinite potentials
        if (beta_u > 0) {
            // Repulsive infinite potential: exp(βu) → ∞, so (1 - exp(βu)) → -∞
            return -g_r;
        } else {
            // Attractive infinite potential: exp(βu) → 0, so (1 - exp(βu)) → 1
            return g_r;
        }
    } else {
        // Normal case with finite potential
        const exp_beta_u = @exp(beta_u);
        return (1.0 - exp_beta_u) * g_r;
    }
}

/// Apply the specified closure relation
pub fn applyClosure(closure_type: ClosureType, g_r: f64, h_r: f64, beta_u: f64) f64 {
    return switch (closure_type) {
        .hypernetted_chain => applyHNC(g_r, h_r, beta_u),
        .percus_yevick => applyPY(g_r, beta_u),
    };
}

// Tests for closure relations
test "PY closure with finite potential" {
    const g_r = 1.0;
    const beta_u = 2.0; // exp(2.0) ≈ 7.39

    const c_r = applyPY(g_r, beta_u);
    const expected = (1.0 - @exp(2.0)) * 1.0; // ≈ -6.39

    try std.testing.expect(@abs(c_r - expected) < 1e-10);
}

test "PY closure with infinite repulsive potential" {
    const g_r = 0.5;
    const beta_u = std.math.inf(f64);

    const c_r = applyPY(g_r, beta_u);
    try std.testing.expectEqual(-g_r, c_r);
}

test "HNC closure with finite potential" {
    const g_r = 2.0;
    const h_r = 1.0;
    const beta_u = 0.5;

    const c_r = applyHNC(g_r, h_r, beta_u);
    const expected = h_r - @log(g_r) - beta_u; // 1.0 - ln(2.0) - 0.5 ≈ -0.193

    try std.testing.expect(@abs(c_r - expected) < 1e-10);
}

test "HNC closure with problematic g_r" {
    const g_r = 0.0; // Problematic for ln(g_r)
    const h_r = -1.0;
    const beta_u = 1.0;

    const c_r = applyHNC(g_r, h_r, beta_u);
    try std.testing.expectEqual(h_r, c_r); // Should fall back to h_r
}
