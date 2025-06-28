const std = @import("std");

/// Convergence acceleration methods
pub const ConvergenceMethod = enum {
    simple_mixing, // Basic linear mixing
    adaptive_mixing, // Adaptive mixing factor
    anderson, // Anderson acceleration
};

/// Convergence control parameters
pub const ConvergenceParams = struct {
    method: ConvergenceMethod = .adaptive_mixing,
    initial_mixing: f64 = 0.1,
    min_mixing: f64 = 0.01,
    max_mixing: f64 = 0.5,
    anderson_depth: usize = 5, // Number of previous iterations to use
    error_increase_factor: f64 = 2.0, // Reduce mixing if error increases by this factor
};

/// Convergence acceleration manager
pub const ConvergenceAccelerator = struct {
    params: ConvergenceParams,
    current_mixing: f64,
    iteration_history: ?[][]f64, // For Anderson acceleration
    error_history: ?[]f64, // Error tracking
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, params: ConvergenceParams) Self {
        return Self{
            .params = params,
            .current_mixing = params.initial_mixing,
            .iteration_history = null,
            .error_history = null,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        if (self.iteration_history) |history| {
            for (history) |h| {
                self.allocator.free(h);
            }
            self.allocator.free(history);
        }
        if (self.error_history) |errors| {
            self.allocator.free(errors);
        }
    }

    /// Initialize Anderson acceleration if needed
    pub fn initAnderson(self: *Self, n_points: usize) !void {
        if (self.params.method != .anderson) return;
        if (self.iteration_history != null) return; // Already initialized

        const depth = self.params.anderson_depth;

        // Allocate history arrays
        const history = try self.allocator.alloc([]f64, depth);
        for (0..depth) |i| {
            history[i] = try self.allocator.alloc(f64, n_points);
        }
        self.iteration_history = history;

        const errors = try self.allocator.alloc(f64, depth);
        self.error_history = errors;
    }

    /// Apply convergence acceleration to values
    pub fn apply(self: *Self, values: []f64, old_values: []f64, iteration: usize) !void {
        switch (self.params.method) {
            .simple_mixing => {
                // Basic linear mixing: h_new = α*h_new + (1-α)*h_old
                const alpha = self.params.initial_mixing;
                for (0..values.len) |i| {
                    const new_value = alpha * values[i] + (1.0 - alpha) * old_values[i];
                    values[i] = if (std.math.isFinite(new_value)) new_value else old_values[i];
                }
            },
            .adaptive_mixing => {
                // Adaptive mixing with current mixing factor
                const alpha = self.current_mixing;
                for (0..values.len) |i| {
                    const new_value = alpha * values[i] + (1.0 - alpha) * old_values[i];
                    values[i] = if (std.math.isFinite(new_value)) new_value else old_values[i];
                }
            },
            .anderson => {
                // Anderson acceleration (simplified implementation)
                if (iteration < self.params.anderson_depth) {
                    // Not enough history yet, use simple mixing
                    const alpha = self.params.initial_mixing;
                    for (0..values.len) |i| {
                        values[i] = alpha * values[i] + (1.0 - alpha) * old_values[i];
                    }

                    // Store in history
                    if (self.iteration_history) |history| {
                        const idx = iteration % self.params.anderson_depth;
                        @memcpy(history[idx], values);
                    }
                } else {
                    // Apply Anderson acceleration
                    try self.andersonAcceleration(values, old_values, iteration);
                }
            },
        }
    }

    /// Update adaptive mixing factor based on error trend
    pub fn updateAdaptiveMixing(self: *Self, current_error: f64, previous_error: f64) void {
        if (self.params.method != .adaptive_mixing) return;

        if (current_error > previous_error * self.params.error_increase_factor) {
            // Error increased significantly, reduce mixing
            self.current_mixing = @max(self.current_mixing * 0.5, self.params.min_mixing);
        } else if (current_error < previous_error * 0.9) {
            // Error decreased, slightly increase mixing
            self.current_mixing = @min(self.current_mixing * 1.1, self.params.max_mixing);
        }
    }

    /// Simplified Anderson acceleration
    fn andersonAcceleration(self: *Self, values: []f64, old_values: []f64, iteration: usize) !void {
        const history = self.iteration_history orelse return;
        const depth = @min(iteration, self.params.anderson_depth);
        const idx = iteration % self.params.anderson_depth;

        // Store current iteration
        @memcpy(history[idx], values);

        // Simple Anderson mixing: average of last few iterations with adaptive weights
        @memset(values, 0.0);

        var total_weight: f64 = 0.0;
        for (0..depth) |i| {
            const weight = 1.0 / @as(f64, @floatFromInt(i + 1)); // Decreasing weights for older iterations
            total_weight += weight;

            for (0..values.len) |j| {
                values[j] += weight * history[(idx + i) % self.params.anderson_depth][j];
            }
        }

        // Normalize
        for (0..values.len) |i| {
            values[i] /= total_weight;
        }

        // Apply some mixing with old value for stability
        const alpha = 0.1;
        for (0..values.len) |i| {
            values[i] = alpha * values[i] + (1.0 - alpha) * old_values[i];
        }
    }
};
