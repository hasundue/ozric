const std = @import("std");

const c = @cImport({
    @cInclude("ceres.h");
});

// Thread-local data management
pub const ThreadLocal = struct {
    init: c.thread_local_init_t,
    deinit: c.thread_local_deinit_t,
};

// Opaque wrapper types
pub const Problem = struct {
    ptr: *anyopaque,

    pub fn create() Problem {
        return Problem{ .ptr = c.ceres_create_problem().? };
    }

    pub fn addResidualBlock(
        self: Problem,
        cost_function: CostFunction,
        parameter_blocks: []*f64,
    ) void {
        // Convert Zig slice of pointers directly to double**
        c.ceres_problem_add_residual_block(
            self.ptr,
            cost_function.ptr,
            @ptrCast(parameter_blocks.ptr),
            @intCast(parameter_blocks.len),
        );
    }

    pub fn solve(self: Problem) bool {
        const result = c.ceres_problem_solve(self.ptr);
        return result == 1;
    }
};

pub const CostFunction = struct {
    ptr: *anyopaque,

    pub fn create(
        cost_fn: c.cost_function_t,
        thread_local: ThreadLocal,
        num_points: c_int,
    ) CostFunction {
        return CostFunction{ .ptr = c.ceres_create_cost_function(
            cost_fn,
            thread_local.init,
            thread_local.deinit,
            num_points,
        ).? };
    }

    pub fn from(ptr: *anyopaque) CostFunction {
        return CostFunction{ .ptr = ptr };
    }
};

test "hello-world optimization" {
    const problem = Problem.create();

    const cost_fn = CostFunction.from(c.ceres_create_hello_world_cost_function().?);

    var x = [_]f64{5.0}; // Initial value: x = 5.0
    var params = [_]*f64{&x[0]};

    // Add the residual block
    problem.addResidualBlock(cost_fn, params[0..]);

    // Solve the optimization problem
    const success = problem.solve();
    try std.testing.expect(success);

    // Check that x converged to ~10.0 (optimal solution)
    const final_x = x[0];
    try std.testing.expect(@abs(final_x - 10.0) < 1e-6);

    std.debug.print("Initial x: 5.0, Final x: {d}\n", .{final_x});
}
