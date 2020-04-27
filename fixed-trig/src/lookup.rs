//! Implementations based on lookup tables.

use fixed::types::I48F16;
use std::convert::TryInto;
use std::include_bytes;

const SIN_TABLE: &[u8] = include_bytes!("tables/sin.table");
const SIN_TABLE_BITS: u8 = SIN_TABLE[0];
const SIN_ANGLE_BITS: u8 = 16;
const QUADRANT_HIGH_MASK: i64 = 1 << (SIN_ANGLE_BITS - 1);
const QUADRANT_LOW_MASK: i64 = 1 << (SIN_ANGLE_BITS - 2);
const SIN_TABLE_OFFSET: u8 = SIN_ANGLE_BITS - 2;
const SIN_INTERP_OFFSET: u8 = SIN_TABLE_OFFSET - SIN_TABLE_BITS;
const SIN_TABLE_MASK: i64 = ((1 << SIN_TABLE_BITS) - 1) << SIN_INTERP_OFFSET;
const SIN_INTERP_MASK: i64 = (1 << SIN_INTERP_OFFSET) - 1;

fn lookup_table(table: &[u8], index: i64) -> i64 {
    let i = (index * 2 + 1) as usize;
    i64::from(u16::from_le_bytes(table[i..(i + 2)].try_into().unwrap()))
}

pub fn sin(angle: I48F16) -> I48F16 {
    let two_pi = I48F16::from_num(2_f64 * std::f64::consts::PI);
    let a = angle.rem_euclid(two_pi);
    let an = (i64::from(a.to_bits()) * (1 << SIN_ANGLE_BITS)) / two_pi.to_bits();
    let is_negative_quadrant: bool = (an & QUADRANT_HIGH_MASK) != 0;
    let is_odd_quadrant: bool = (an & QUADRANT_LOW_MASK) == 0;
    let mut index = (an & SIN_TABLE_MASK) >> SIN_INTERP_OFFSET;
    if !is_odd_quadrant {
        index = (1 << SIN_TABLE_BITS) - 1 - index;
    }
    let x1 = lookup_table(SIN_TABLE, index);
    let x2 = lookup_table(SIN_TABLE, index + 1);
    let interp = an & SIN_INTERP_MASK;
    let approx = ((x2 - x1) * interp) >> SIN_INTERP_OFFSET;
    let mut sine = if is_odd_quadrant {
        x1 + approx
    } else {
        x2 - approx
    };
    if is_negative_quadrant {
        sine = -sine;
    }
    I48F16::from_bits(sine)
}

pub fn cos(val: I48F16) -> I48F16 {
    let cos_offset = I48F16::from_num((5_f64 * std::f64::consts::PI) / 2_f64);
    sin(val + cos_offset)
}

const ASIN_TABLE: &[u8] = include_bytes!("tables/asin.table");
const ASIN_TABLE_BITS: u8 = ASIN_TABLE[0];
const ASIN_VAL_BITS: u8 = 16;
const ASIN_TABLE_OFFSET: u8 = ASIN_VAL_BITS;
const ASIN_INTERP_OFFSET: u8 = ASIN_TABLE_OFFSET - ASIN_TABLE_BITS;
const ASIN_TABLE_MASK: i64 = ((1 << ASIN_TABLE_BITS) - 1) << ASIN_INTERP_OFFSET;
const ASIN_INTERP_MASK: i64 = (1 << ASIN_INTERP_OFFSET) - 1;

pub fn asin(mut val: I48F16) -> I48F16 {
    let pi_2 = I48F16::from_num(std::f64::consts::FRAC_PI_2);
    let is_neg = val < 0;
    if is_neg {
        val = -val;
    }
    if val > 1 {
        panic!("invalid asin value");
    } else if val == 1 {
        return if is_neg { -pi_2 } else { pi_2 };
    }
    let val_bits = val.to_bits();
    let index = (val_bits & ASIN_TABLE_MASK) >> ASIN_INTERP_OFFSET;
    let x1 = lookup_table(ASIN_TABLE, index);
    let x2 = lookup_table(ASIN_TABLE, index + 1);
    let interp = val_bits & ASIN_INTERP_MASK;
    let approx = ((x2 - x1) * interp) >> ASIN_INTERP_OFFSET;
    let mut asine = x1 + approx;
    if is_neg {
        asine = -asine;
    }
    I48F16::from_bits(asine) * pi_2
}

pub fn acos(val: I48F16) -> I48F16 {
    asin(-val) + I48F16::from_num(std::f64::consts::FRAC_PI_2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin() {
        for i in -100..100 {
            let fx = f64::from(i) * 0.1_f64;
            let x: I48F16 = I48F16::from_num(fx);
            let sin_x: f64 = sin(x).to_num();
            let sin_fx = fx.sin();
            let dsx = (sin_x - sin_fx).abs();
            if dsx > 0.001 {
                panic!("mismatch sin({}): {} {}", x, sin_x, sin_fx);
            }
        }
    }

    #[test]
    fn test_sin_comprehensive() {
        for i in 0..(1 << 20) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            let sin_x: f64 = sin(x).to_num();
            let sin_fx = fx.sin();
            let dsx = (sin_x - sin_fx).abs();
            if dsx > 0.001 {
                panic!("mismatch sin({}): {} {}", x, sin_x, sin_fx);
            }
        }
    }

    #[test]
    fn test_cos() {
        for i in -100..100 {
            let fx = f64::from(i) * 0.1_f64;
            let x: I48F16 = I48F16::from_num(fx);
            let cos_x: f64 = cos(x).to_num();
            let cos_fx = fx.cos();
            let dsx = (cos_x - cos_fx).abs();
            if dsx > 0.001 {
                panic!("mismatch cos({}): {} {}", x, cos_x, cos_fx);
            }
        }
    }

    #[test]
    fn test_cos_comprehensive() {
        for i in 0..(1 << 20) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            let cos_x: f64 = cos(x).to_num();
            let cos_fx = fx.cos();
            let dsx = (cos_x - cos_fx).abs();
            if dsx > 0.001 {
                panic!("mismatch cos({}): {} {}", x, cos_x, cos_fx);
            }
        }
    }

    #[test]
    fn test_asin_comprehensive() {
        let one = 1 << 16;
        for i in -one..(one + 1) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            let asin_x: f64 = asin(x).to_num();
            let asin_fx = fx.asin();
            let dsx = (asin_x - asin_fx).abs();
            // does not do well at values very close to 1 or -1
            if dsx > 0.001 && ((fx > -0.999 && fx < 0.999) || dsx > 0.02) {
                panic!("mismatch asin({}): {} {}", x, asin_x, asin_fx);
            }
        }
    }

    #[test]
    fn test_acos_comprehensive() {
        let one = 1 << 16;
        for i in -one..(one + 1) {
            let x = I48F16::from_bits(i);
            let fx: f64 = x.to_num();
            let acos_x: f64 = acos(x).to_num();
            let acos_fx = fx.acos();
            let dsx = (acos_x - acos_fx).abs();
            // does not do well at values very close to 1 or -1
            if dsx > 0.001 && ((fx > -0.999 && fx < 0.999) || dsx > 0.02) {
                panic!("mismatch acos({}): {} {}", x, acos_x, acos_fx);
            }
        }
    }
}
