//! Implementations of special functions based on the CORDIC algorithm.

mod cordic_number;

pub use cordic_number::CordicNumber;
use fixed::types::U0F64;
use std::convert::TryInto;

const ATAN_TABLE: &[u8] = include_bytes!("tables/cordic_atan.table");
const EXP_MINUS_ONE_TABLE: &[u8] = include_bytes!("tables/cordic_exp_minus_one.table");

fn lookup_table(table: &[u8], index: u8) -> U0F64 {
    let i = index as usize * 8;
    U0F64::from_bits(u64::from_le_bytes(table[i..(i + 8)].try_into().unwrap()))
}

// See cordit1 from http://www.voidware.com/cordic.htm
fn cordic_circular<T: CordicNumber>(mut x: T, mut y: T, mut z: T, vecmode: T) -> (T, T, T) {
    let _0 = T::zero();
    let _2 = T::one() + T::one();

    for i in 0..T::num_fract_bits() {
        if vecmode >= _0 && y < vecmode || vecmode < _0 && z >= _0 {
            let x1 = x - (y >> i);
            y = y + (x >> i);
            x = x1;
            z = z - T::from_u0f64(lookup_table(ATAN_TABLE, i));
        } else {
            let x1 = x + (y >> i);
            y = y - (x >> i);
            x = x1;
            z = z + T::from_u0f64(lookup_table(ATAN_TABLE, i));
        }
    }

    (x, y, z)
}

fn gain_cordic<T: CordicNumber>() -> T {
    cordic_circular(T::one(), T::zero(), T::zero(), -T::one()).0
}

/// Compute simultaneously the sinus and cosine of the given fixed-point number.
pub fn sin_cos<T: CordicNumber>(mut angle: T) -> (T, T) {
    let mut negative = false;

    while angle > T::frac_pi_2() {
        angle -= T::pi();
        negative = !negative;
    }

    while angle < -T::frac_pi_2() {
        angle += T::pi();
        negative = !negative;
    }

    let inv_gain = T::one() / gain_cordic(); // FIXME: precompute this.
    let res = cordic_circular(inv_gain, T::zero(), angle, -T::one());

    if negative {
        (-res.1, -res.0)
    } else {
        (res.1, res.0)
    }
}

/// Compute the sinus of the given fixed-point number.
pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

/// Compute the cosinus of the given fixed-point number.
pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

/// Compute the tangent of the given fixed-point number.
pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (sin, cos) = sin_cos(angle);
    sin / cos
}

/// Compute the arc-sinus of the given fixed-point number.
pub fn asin<T: CordicNumber>(mut val: T) -> T {
    // For asin, we use a double-rotation approach to reduce errors.
    // NOTE: see https://stackoverflow.com/questions/25976656/cordic-arcsine-implementation-fails
    // for details about the innacuracy of CORDIC for asin.

    let mut theta = T::zero();
    let mut z = (T::one(), T::zero());
    let niter = T::num_fract_bits();

    for j in 0..niter {
        let sign_x = if z.0 < T::zero() { -T::one() } else { T::one() };
        let sigma = if z.1 <= val { sign_x } else { -sign_x };
        let rotate = |(x, y)| (x - ((y >> j) * sigma), y + ((x >> j) * sigma));
        z = rotate(rotate(z));

        let angle = T::from_u0f64(lookup_table(ATAN_TABLE, j));
        theta = theta + ((angle + angle) * sigma);
        val = val + (val >> (j + j));
    }

    theta
}

/// Compute the arc-cosine of the given fixed-point number.
pub fn acos<T: CordicNumber>(val: T) -> T {
    T::frac_pi_2() - asin(val)
}

/// Compute the arc-tangent of the given fixed-point number.
pub fn atan<T: CordicNumber>(val: T) -> T {
    cordic_circular(T::one(), val, T::zero(), T::zero()).2
}

/// Compute the arc-tangent of `y/x` with quadrant correction.
pub fn atan2<T: CordicNumber>(y: T, x: T) -> T {
    if x == T::zero() {
        if y < T::zero() {
            return -T::frac_pi_2();
        } else {
            return T::frac_pi_2();
        }
    }

    if y == T::zero() {
        if x >= T::zero() {
            return T::zero();
        } else {
            return T::pi();
        }
    }

    match (x < T::zero(), y < T::zero()) {
        (false, false) => atan(y / x),
        (false, true) => -atan(-y / x),
        (true, false) => T::pi() - atan(y / -x),
        (true, true) => atan(y / x) - T::pi(),
    }
}

/// Compute the exponential root of the given fixed-point number.
pub fn exp<T: CordicNumber>(x: T) -> T {
    assert!(
        T::num_fract_bits() <= 128,
        "Exp is not supported for more than 128 decimals."
    );
    let _0 = T::zero();
    let _1 = T::one();
    let _3 = T::one() + T::one() + T::one();
    let mut int_part = x.floor();
    let mut dec_part = x - int_part;
    let mut poweroftwo = T::half();
    let mut w = [false; 128];

    for i in 0..T::num_fract_bits() {
        if poweroftwo < dec_part {
            w[i as usize] = true;
            dec_part -= poweroftwo;
        }

        poweroftwo = poweroftwo >> 1;
    }

    let mut fx = _1;

    for i in 0..T::num_fract_bits() {
        if w[i as usize] {
            let ai = T::from_u0f64(lookup_table(EXP_MINUS_ONE_TABLE, i)) + T::one();
            fx = fx * ai;
        }
    }

    let f4 = _1 + (dec_part >> 2);
    let f3 = _1 + (dec_part / _3) * f4;
    let f2 = _1 + (dec_part >> 1) * f3;
    let f1 = _1 + dec_part * f2;
    fx = fx * f1;

    if int_part < _0 {
        while int_part != _0 {
            fx = fx / T::e();
            int_part += _1;
        }
    } else {
        while int_part != _0 {
            fx = fx * T::e();
            int_part -= _1;
        }
    }

    fx
}

/// Compute the square root of the given fixed-point number.
pub fn sqrt<T: CordicNumber>(x: T) -> T {
    if x == T::zero() || x == T::one() {
        return x;
    }

    let mut pow2 = T::one();
    let mut result;

    if x < T::one() {
        while x <= pow2 * pow2 {
            pow2 = pow2 >> 1;
        }

        result = pow2;
    } else {
        // x >= T::one()
        while pow2 * pow2 <= x {
            pow2 = pow2 << 1;
        }

        result = pow2 >> 1;
    }

    for _ in 0..T::num_bits() {
        pow2 = pow2 >> 1;
        let next_result = result + pow2;
        if next_result * next_result <= x {
            result = next_result;
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use fixed::types::I48F16;

    fn assert_approx_eq<T: std::fmt::Display>(
        input: T,
        computed: f64,
        expected: f64,
        max_err: f64,
    ) {
        let err = (computed - expected).abs();
        if err > max_err {
            panic!(
                "mismatch for input {}: computed {}, expected {}",
                input, computed, expected
            );
        }
    }

    macro_rules! test_trig(
        ($test: ident, $test_comprehensive: ident, $trigf: ident, $max_err: expr) => {
            #[test]
            fn $test() {
                for i in -100..100 {
                    let fx = f64::from(i) * 0.1_f64;
                    let x: I48F16 = I48F16::from_num(fx);
                    assert_approx_eq(x, $trigf(x).to_num(), fx.$trigf(), $max_err);
                }
            }

            #[test]
            fn $test_comprehensive() {
                for i in 0..(1 << 20) {
                    let x = I48F16::from_bits(i);
                    let fx: f64 = x.to_num();
                    assert_approx_eq(x, $trigf(x).to_num(), fx.$trigf(), $max_err);

                    // Test negative numbers too.
                    let x = -I48F16::from_bits(i);
                    let fx: f64 = x.to_num();
                    assert_approx_eq(x, $trigf(x).to_num(), fx.$trigf(), $max_err);
                }
            }
        }
    );

    test_trig!(test_sin, test_sin_comprehensive, sin, 0.001);
    test_trig!(test_cos, test_cos_comprehensive, cos, 0.001);
    test_trig!(test_atan, test_atan_comprehensive, atan, 0.001);

    #[test]
    fn test_asin() {
        for i in 0..(1 << 17) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, asin(x).to_num(), fx.asin(), 0.01);

            // Test negative numbers too.
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, asin(x).to_num(), fx.asin(), 0.01);
        }
    }

    #[test]
    fn test_acos() {
        for i in 0..(1 << 17) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, acos(x).to_num(), fx.acos(), 0.01);

            // Test negative numbers too.
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, acos(x).to_num(), fx.acos(), 0.01);
        }
    }

    #[test]
    fn test_sqrt() {
        for i in 0..(1 << 20) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, sqrt(x).to_num(), fx.sqrt(), 0.01);
        }
    }

    #[test]
    fn test_exp() {
        for i in 0..(1 << 18) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            assert_approx_eq(x, exp(x).to_num(), fx.exp(), 0.01);
        }
    }
}
