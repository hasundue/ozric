const std = @import("std");

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
