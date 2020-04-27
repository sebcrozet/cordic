use fixed::types::extra::{
    IsLessOrEqual, True, Unsigned, U13, U14, U16, U29, U30, U32, U5, U6, U61, U62, U64, U8,
};
use fixed::types::U0F64;
use fixed::{FixedI16, FixedI32, FixedI64, FixedI8};
use std::ops::{Add, AddAssign, BitXor, Div, Mul, MulAssign, Neg, Shl, Shr, Sub, SubAssign};

pub trait CordicNumber:
    Copy
    + PartialOrd
    + AddAssign
    + SubAssign
    + Div<Output = Self>
    + Mul<Output = Self> // FIXME: remove that
    + Neg<Output = Self>
    + Sub<Output = Self>
    + Add<Output = Self>
    + Shr<u8, Output = Self>
    + Shl<u8, Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn frac_pi_2() -> Self;
    fn pi() -> Self;
    fn from_u0f64(val: U0F64) -> Self;
    fn num_fract_bits() -> u8;
}

// The IsLessOrEqual constraints are for (in order):
// - The Fixed type
// - The FRAC_PI_2 constant.
// - The PI constant.
impl<Fract> CordicNumber for FixedI8<Fract>
where
    Fract: Unsigned
        + IsLessOrEqual<U8, Output = True>
        + IsLessOrEqual<U6, Output = True>
        + IsLessOrEqual<U5, Output = True>,
{
    #[inline(always)]
    fn zero() -> Self {
        Self::from_bits(0)
    }

    #[inline(always)]
    fn one() -> Self {
        Self::from_num(1.0)
    }

    #[inline(always)]
    fn frac_pi_2() -> Self {
        Self::FRAC_PI_2
    }

    #[inline(always)]
    fn pi() -> Self {
        Self::PI
    }

    #[inline(always)]
    fn from_u0f64(val: U0F64) -> Self {
        val.to_num()
    }

    #[inline(always)]
    fn num_fract_bits() -> u8 {
        Fract::to_u8()
    }
}

impl<Fract> CordicNumber for FixedI32<Fract>
where
    Fract: Unsigned
        + IsLessOrEqual<U32, Output = True>
        + IsLessOrEqual<U30, Output = True>
        + IsLessOrEqual<U29, Output = True>,
{
    #[inline(always)]
    fn zero() -> Self {
        Self::from_bits(0)
    }

    #[inline(always)]
    fn one() -> Self {
        Self::from_num(1.0)
    }

    #[inline(always)]
    fn frac_pi_2() -> Self {
        Self::FRAC_PI_2
    }

    #[inline(always)]
    fn pi() -> Self {
        Self::PI
    }

    #[inline(always)]
    fn from_u0f64(val: U0F64) -> Self {
        val.to_num()
    }

    #[inline(always)]
    fn num_fract_bits() -> u8 {
        Fract::to_u8()
    }
}

impl<Fract> CordicNumber for FixedI16<Fract>
where
    Fract: Unsigned
        + IsLessOrEqual<U16, Output = True>
        + IsLessOrEqual<U14, Output = True>
        + IsLessOrEqual<U13, Output = True>,
{
    #[inline(always)]
    fn zero() -> Self {
        Self::from_bits(0)
    }

    #[inline(always)]
    fn one() -> Self {
        Self::from_num(1.0)
    }

    #[inline(always)]
    fn frac_pi_2() -> Self {
        Self::FRAC_PI_2
    }

    #[inline(always)]
    fn pi() -> Self {
        Self::PI
    }

    #[inline(always)]
    fn from_u0f64(val: U0F64) -> Self {
        val.to_num()
    }

    #[inline(always)]
    fn num_fract_bits() -> u8 {
        Fract::to_u8()
    }
}

impl<Fract> CordicNumber for FixedI64<Fract>
where
    Fract: Unsigned
        + IsLessOrEqual<U64, Output = True>
        + IsLessOrEqual<U62, Output = True>
        + IsLessOrEqual<U61, Output = True>,
{
    #[inline(always)]
    fn zero() -> Self {
        Self::from_bits(0)
    }

    #[inline(always)]
    fn one() -> Self {
        Self::from_num(1.0)
    }

    #[inline(always)]
    fn frac_pi_2() -> Self {
        Self::FRAC_PI_2
    }

    #[inline(always)]
    fn pi() -> Self {
        Self::PI
    }

    #[inline(always)]
    fn from_u0f64(val: U0F64) -> Self {
        val.to_num()
    }

    #[inline(always)]
    fn num_fract_bits() -> u8 {
        Fract::to_u8()
    }
}