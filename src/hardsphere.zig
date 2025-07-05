const std = @import("std");
const math = std.math;
const t = std.testing;

pub fn weightFunctionZeroth(
    sigma: f64,
    r: f64,
) f64 {
    if (r > sigma) return 0.0;

    const sigma3 = math.pow(f64, sigma, 3);
    const pi = math.pi;

    return 3 / (4 * pi * sigma3);
}

pub fn weightFunctionFirst(
    sigma: f64,
    r: f64,
) f64 {
    const x = r / sigma;
    const x2 = math.pow(f64, x, 2);

    if (x <= 1.0) {
        const a_0 = 0.90724;
        const a_1 = -1.23717;
        const a_2 = 0.21616;

        return a_0 + a_1 * x + a_2 * x2;
    } else {
        const x3 = math.pow(f64, x, 3);

        const c = -0.10244;
        const b_0 = 35.134;
        const alpha = 4.934;
        const b_1 = -98.684;
        const beta_1 = 3.5621;
        const b_2 = 92.693;
        const beta_2 = 12.0;
        const b_3 = -29.257;

        const term1 = c * @exp(-beta_1 * (x - 1.0)) * @sin(alpha * (x - 1.0));
        const term2 = @exp(-beta_2 * (x - 1.0)) * (b_0 + b_1 * x + b_2 * x2 + b_3 * x3);
        return term1 + term2;
    }
}

test weightFunctionFirst {
    const sigma = 1.0;

    const w1 = struct {
        pub fn call(r: f64) f64 {
            return weightFunctionFirst(sigma, r);
        }
    }.call;

    try t.expectApproxEqAbs(0.90724, w1(0), 1e-5);

    // Should be positive deep inside the sphere
    try t.expect(w1(0.5) > 0);

    // Should be negative at contact distance
    try t.expect(w1(1.0) < 0);

    // Test continuity at contact distance
    const eps = 1e-6;
    try t.expectApproxEqAbs(w1(sigma - eps), w1(sigma + eps), 1e-3);

    // Test oscillating structure of the tail
    try t.expect(w1(1.5) < 0);
    try t.expect(w1(2.0) > 0);
    try t.expect(w1(2.5) < 0);
}

pub fn weightFunctionSecond(
    sigma: f64,
    r: f64,
) f64 {
    if (r >= sigma) return 0.0;

    const x = r / sigma;
    const x2 = x * x;
    const pi = math.pi;

    return (15 / 8 / pi / sigma) * (1 - 3 * x + 3 * x2);
}
