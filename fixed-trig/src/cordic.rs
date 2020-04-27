//! Implementations based on the CORDIC algorithm.
pub use super::cordic_number::CordicNumber;
use fixed::types::U0F64;
use std::convert::TryInto;

const ATAN_TABLE: &[u8] = include_bytes!("tables/cordic_atan.table");

fn lookup_table(table: &[u8], index: u8) -> U0F64 {
    let i = index as usize * 8;
    U0F64::from_bits(u64::from_le_bytes(table[i..(i + 8)].try_into().unwrap()))
}

pub fn cordic_circular<T: CordicNumber>(mut x: T, mut y: T, mut z: T, vecmode: T) -> (T, T, T) {
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

pub fn gain_cordic<T: CordicNumber>() -> T {
    cordic_circular(T::one(), T::zero(), T::zero(), -T::one()).0
}

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

pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (sin, cos) = sin_cos(angle);
    sin / cos
}

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

pub fn acos<T: CordicNumber>(val: T) -> T {
    T::frac_pi_2() - asin(val)
}

pub fn atan<T: CordicNumber>(val: T) -> T {
    cordic_circular(T::one(), val, T::zero(), T::zero()).2
}

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

// FIXME: this should be contributed to the fixed-sqrt crate instead of remaining here.
pub fn sqrt<T: CordicNumber>(x: T, niters: usize) -> T {
    if x == T::zero() || x == T::one() {
        return x;
    }

    // FIXME: optimize with bitshifts
    let mut pow2 = T::one();
    let mut result = T::zero();

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

    for _ in 0..niters {
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
            assert_approx_eq(x, sqrt(x, 40).to_num(), fx.sqrt(), 0.01);
        }
    }
}
