//! Rust version by Rik Huijzer and Jose Storopoli

use crate::dpq::*;
use crate::libc::*;
use crate::nmath::*;
use crate::rmath::*;

/// Compute the quantile function for the normal distribution.
///
/// The algorithm AS 241 of Wichura is used,
///     and has been improved for the very extreme tail (and log_p=TRUE)
///
/// ## REFERENCE
///
/// Wichura, M.J. (1988).
/// Algorithm AS 241: The Percentage Points of the Normal Distribution.
/// Applied Statistics, 37, 477-484.
pub fn qnorm5(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    let mut r: f64;
    let mut val: f64;

    if p.is_nan() || mu.is_nan() || sigma.is_nan() {
        return p + mu + sigma;
    }

    if let Some(x) = r_q_p01_boundaries(p, ML_NEGINF, ML_POSINF, lower_tail, log_p) {
        return x;
    }

    if sigma < 0.0 {
        ml_warn_return_nan();
    }
    if sigma == 0.0 {
        return mu;
    }

    let p_: f64 = r_dt_qiv(p, lower_tail, log_p); // real lower_tail prob. p
    let q: f64 = p_ - 0.5;

    // - use AS 241 --- */
    // double ppnd16_(double *p, long *ifault)*/
    //      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
    //
    //      Produces the normal deviate Z corresponding to a given lower
    //      tail area of P; Z is accurate to about 1 part in 10**16.
    //
    //      (original fortran code used PARAMETER(..) for the coefficients
    //       and provided hash codes for checking them...)
    //
    if q.abs() <= 0.425 {
        // |p~ - 0.5| <= .425  <==> 0.075 <= p~ <= 0.925
        r = 0.180625 - q * q; // = .425^2 - q^2  >= 0
        val = q
            * (((((((r * 2_509.080_928_730_122_7 + 33_430.575_583_588_13) * r
                + 67_265.770_927_008_7)
                * r
                + 45_921.953_931_549_87)
                * r
                + 13_731.693_765_509_46)
                * r
                + 1_971.590_950_306_551_3)
                * r
                + 133.141_667_891_784_38)
                * r
                + 3.387_132_872_796_366_5)
            / (((((((r * 5_226.495_278_852_854 + 28_729.085_735_721_943) * r
                + 39_307.895_800_092_71)
                * r
                + 21_213.794_301_586_597)
                * r
                + 5_394.196_021_424_751)
                * r
                + 687.187_007_492_057_9)
                * r
                + 42.313_330_701_600_91)
                * r
                + 1.);
    } else {
        /* closer than 0.075 from {0,1} boundary :
         *  r := log(p~);  p~ = min(p, 1-p) < 0.075 :  */
        if log_p && ((lower_tail && q <= 0.0) || (!lower_tail && q > 0.0)) {
            r = p;
        } else {
            let val = if q > 0.0 {
                r_dt_civ(p, lower_tail, log_p) /* 1-p */
            } else {
                p_ /* = R_DT_Iv(p) ^=  p */
            };
            r = log(val);
        }
        // r = sqrt( - log(min(p,1-p)) )  <==>  min(p, 1-p) = exp( - r^2 ) :
        r = sqrt(-r);
        if r <= 5. {
            // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11
            r += -1.6;
            val = (((((((r * 7.745_450_142_783_414e-4 + 0.022_723_844_989_269_184) * r
                + 0.241_780_725_177_450_6)
                * r
                + 1.270_458_252_452_368_4)
                * r
                + 3.647_848_324_763_204_5)
                * r
                + 5.769_497_221_460_691)
                * r
                + 4.630_337_846_156_546)
                * r
                + 1.423_437_110_749_683_5)
                / (((((((r * 1.050_750_071_644_416_9e-9 + 5.475_938_084_995_345e-4) * r
                    + 0.015_198_666_563_616_457)
                    * r
                    + 0.148_103_976_427_480_08)
                    * r
                    + 0.689_767_334_985_1)
                    * r
                    + 1.676_384_830_183_803_8)
                    * r
                    + 2.053_191_626_637_759)
                    * r
                    + 1.);
        } else if r >= 816.0 {
            // p is *extremly* close to 0 or 1 - only possibly when log_p =TRUE
            // Using the asymptotical formula -- is *not* optimal but uniformly better than branch below
            val = r * M_SQRT2;
        } else {
            // p is very close to  0 or 1:  r > 5 <==> min(p,1-p) < exp(-25) = 1.3888..e-11
            // Wichura, p.478: minimax rational approx R_3(t) is for 5 <= t <= 27  (t :== r)
            r += -5.;
            val = (((((((r * 2.010_334_399_292_288_1e-7 + 2.711_555_568_743_487_6e-5) * r
                + 0.001_242_660_947_388_078_4)
                * r
                + 0.026_532_189_526_576_124)
                * r
                + 0.296_560_571_828_504_87)
                * r
                + 1.784_826_539_917_291_3)
                * r
                + 5.463_784_911_164_114)
                * r
                + 6.657_904_643_501_103)
                / (((((((r * 2.044_263_103_389_939_7e-15 + 1.421_511_758_316_446e-7) * r
                    + 1.846_318_317_510_054_8e-5)
                    * r
                    + 7.868_691_311_456_133e-4)
                    * r
                    + 0.014_875_361_290_850_615)
                    * r
                    + 0.136_929_880_922_735_8)
                    * r
                    + 0.599_832_206_555_888)
                    * r
                    + 1.);
        }

        if q < 0.0 {
            val = -val;
            /* return (q >= 0.)? r : -r ;*/
        }
    }
    mu + sigma * val
}
