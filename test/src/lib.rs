#[cfg(test)]
mod test_math {
    use approx::abs_diff_eq;
    use rmathlib::*;

    mod pbeta;

    mod c {
        extern "C" {
            pub fn Rf_bd0(x: f64, np: f64) -> f64;
            pub fn Rf_chebyshev_eval(x: f64, a: *mut f64, n: i32) -> f64;
            pub fn Rf_chebyshev_init(dos: *mut f64, nos: i32, eta: f64) -> i32;
            pub fn Rf_ebd0(x: f64, M: f64, yh: *mut f64, yl: *mut f64);
            pub fn Rf_i1mach(i: i32) -> i32;
            pub fn Rf_lgammacor(x: f64) -> f64;
            pub fn Rf_stirlerr(n: f64) -> f64;
            pub fn cospi(x: f64) -> f64;
            pub fn dgamma(x: f64, shape: f64, scale: f64, give_log: bool) -> f64;
            pub fn dnorm4(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64;
            pub fn dpois(x: f64, lambda: f64, give_log: bool) -> f64;
            pub fn gammafn(x: f64) -> f64;
            pub fn lbeta(a: f64, b: f64) -> f64;
            pub fn log1pmx(x: f64) -> f64;
            pub fn lgammafn(x: f64) -> f64;
            pub fn lgammafn_sign(x: f64, sgn: Option<&mut i32>) -> f64;
            pub fn pgamma(x: f64, alph: f64, scale: f64, lower_tail: i32, log_p: i32) -> f64;
            pub fn pnorm5(x: f64, mu: f64, sigma: f64, lower_tail: i32, log_p: i32) -> f64;
            #[allow(dead_code)]
            pub fn pt(x: f64, n: f64, lower_tail: i32, log_p: i32) -> f64;
            pub fn qnorm5(p: f64, mu: f64, sigma: f64, lower_tail: i32, log_p: i32) -> f64;
            pub fn sinpi(x: f64) -> f64;
            pub fn tanpi(x: f64) -> f64;
        }
    }

    #[test]
    fn test_bd0() {
        // if: 9.0-8.9 < 0.1*(9.0+8.9) = 0.1 < 1.79
        assert_eq!(bd0(9.0, 8.9), unsafe { c::Rf_bd0(9.0, 8.9) });
        // else: ((5.0-1.0).abs !< 0.1*(5.0+1.0)) = 4.0 !< 0.6
        assert_eq!(bd0(5.0, 1.0), unsafe { c::Rf_bd0(5.0, 1.0) });
    }

    #[test]
    fn test_chebyshev_eval() {
        assert_eq!(chebyshev_eval(0.6, &[1.0, 2.0, 3.0], 2), unsafe {
            c::Rf_chebyshev_eval(0.6, [1.0, 2.0, 3.0].as_mut_ptr(), 2)
        });
        assert!(chebyshev_eval(0.6, &[1.0, 2.0, 3.0], 0).is_nan());
        assert_eq!(chebyshev_eval(0.6, &[1.0, 2.0, 3.0], 2), unsafe {
            c::Rf_chebyshev_eval(0.6, [1.0, 2.0, 3.0].as_mut_ptr(), 2)
        });
    }

    #[test]
    fn test_chebyshev_init() {
        assert_eq!(chebyshev_init(&[1.0, 2.0, 3.0], 3, 0.5), unsafe {
            c::Rf_chebyshev_init([1.0, 2.0, 3.0].as_mut_ptr(), 3, 0.5)
        });
        assert_eq!(chebyshev_init(&[], 0, 0.5), unsafe {
            c::Rf_chebyshev_init([].as_mut_ptr(), 0, 0.5)
        });
        assert_eq!(chebyshev_init(&[1.0, 2.0, 3.0], 3, -0.5), unsafe {
            c::Rf_chebyshev_init([1.0, 2.0, 3.0].as_mut_ptr(), 3, -0.5)
        });
    }

    fn test_ebd0_helper(x: f64, m: f64) {
        let mut c_yh: f64 = f64::NAN;
        let mut c_yl: f64 = f64::NAN;
        let (yh, yl) = ebd0(x, m);
        unsafe { c::Rf_ebd0(x, m, &mut c_yh, &mut c_yl) };
        assert_eq!(yh, c_yh, "yh with x={x:?}, m={m:?}", x = x, m = m);
        assert_eq!(yl, c_yl, "yl with x={x:?}, m={m:?}", x = x, m = m);
    }

    #[test]
    fn test_ebd0() {
        test_ebd0_helper(1.0, 1.0);
        test_ebd0_helper(1.0, 1.0);
        test_ebd0_helper(0.0, 1.0);
        test_ebd0_helper(1.0, 0.0);
        test_ebd0_helper(3.0, 0.5);
        test_ebd0_helper(10.2, 5.45);
    }

    #[test]
    fn test_i1mach() {
        assert_eq!(i1mach(1), unsafe { c::Rf_i1mach(1) });
        assert_eq!(i1mach(2), unsafe { c::Rf_i1mach(2) });
        assert_eq!(i1mach(3), unsafe { c::Rf_i1mach(3) });
        assert_eq!(i1mach(4), unsafe { c::Rf_i1mach(4) });
        assert_eq!(i1mach(5), unsafe { c::Rf_i1mach(5) });
        assert_eq!(i1mach(6), unsafe { c::Rf_i1mach(6) });
        assert_eq!(i1mach(7), unsafe { c::Rf_i1mach(7) });
        assert_eq!(i1mach(8), unsafe { c::Rf_i1mach(8) });
        assert_eq!(i1mach(9), unsafe { c::Rf_i1mach(9) });
        assert_eq!(i1mach(10), unsafe { c::Rf_i1mach(10) });
        assert_eq!(i1mach(11), unsafe { c::Rf_i1mach(11) });
        assert_eq!(i1mach(12), unsafe { c::Rf_i1mach(12) });
        assert_eq!(i1mach(13), unsafe { c::Rf_i1mach(13) });
        assert_eq!(i1mach(14), unsafe { c::Rf_i1mach(14) });
        assert_eq!(i1mach(15), unsafe { c::Rf_i1mach(15) });
        assert_eq!(i1mach(16), unsafe { c::Rf_i1mach(16) });
    }

    #[test]
    fn test_lgammacor() {
        assert!(lgammacor(-1.0).is_nan());
        assert!(lgammacor(0.0).is_nan());
        assert_eq!(lgammacor(10.0), unsafe { c::Rf_lgammacor(10.0) });
        assert_eq!(lgammacor(20.0), unsafe { c::Rf_lgammacor(20.0) });
        assert_eq!(lgammacor(3.745194030963158e306 + 1.0), unsafe {
            c::Rf_lgammacor(3.745194030963158e306 + 1.0)
        });
        assert_eq!(lgammacor(94906265.62425156 + 1.0), unsafe {
            c::Rf_lgammacor(94906265.62425156 + 1.0)
        });
    }

    #[test]
    fn test_stirlerr() {
        assert_eq!(stirlerr(1.0), unsafe { c::Rf_stirlerr(1.0) });
        assert_eq!(stirlerr(2.0), unsafe { c::Rf_stirlerr(2.0) });
        assert_eq!(stirlerr(25.0), unsafe { c::Rf_stirlerr(25.0) });
        assert_eq!(stirlerr(50.0), unsafe { c::Rf_stirlerr(50.0) });
    }

    #[test]
    fn test_cospi() {
        assert_eq!(cospi(0.0), unsafe { c::cospi(0.0) });
        assert_eq!(cospi(0.234), unsafe { c::cospi(0.234) });
    }

    #[test]
    fn test_dgamma() {
        assert_eq!(dgamma(0.0, 0.0, 1.0, false), unsafe {
            c::dgamma(0.0, 0.0, 1.0, false)
        });
        assert!(
            dgamma(0.0, -1.0, 1.0, false).is_nan()
                && unsafe { c::dgamma(0.0, -1.0, 1.0, false).is_nan() }
        );
        assert_eq!(dgamma(0.12, 0.34, 0.56, true), unsafe {
            c::dgamma(0.12, 0.34, 0.56, true)
        });
        assert_eq!(dgamma(0.12, 0.34, 0.56, false), unsafe {
            c::dgamma(0.12, 0.34, 0.56, false)
        });
        assert_eq!(dgamma(0.12, 1.34, 0.56, false), unsafe {
            c::dgamma(0.12, 1.34, 0.56, false)
        });
    }

    #[test]
    fn test_dnorm() {
        assert_eq!(dnorm(0.0, 0.0, 1.0, false), unsafe {
            c::dnorm4(0.0, 0.0, 1.0, false)
        });
        assert_eq!(dnorm(0.0, 0.0, 1.0, true), unsafe {
            c::dnorm4(0.0, 0.0, 1.0, true)
        });
        assert_eq!(dnorm(1.0, 0.0, 1.0, false), unsafe {
            c::dnorm4(1.0, 0.0, 1.0, false)
        });
        assert_eq!(dnorm(1.0, 0.0, 1.0, true), unsafe {
            c::dnorm4(1.0, 0.0, 1.0, true)
        });
        assert_eq!(dnorm(-1.0, 0.0, 1.0, false), unsafe {
            c::dnorm4(-1.0, 0.0, 1.0, false)
        });
        assert_eq!(dnorm(-1.0, 0.0, 1.0, true), unsafe {
            c::dnorm4(-1.0, 0.0, 1.0, true)
        });
    }

    #[test]
    fn test_dpois() {
        assert!(dpois(-1.0, -1.0, false).is_nan());
        assert!(dpois(1.0, -1.0, false).is_nan());
        assert_eq!(dpois(1.0, 1.0, false), unsafe { c::dpois(1.0, 1.0, false) });
        assert_eq!(dpois(1.0, 1.0, true), unsafe { c::dpois(1.0, 1.0, true) });
    }

    #[test]
    fn test_gammafn() {
        assert!(gammafn(-1.0).is_nan());
        assert!(gammafn(0.0).is_nan());
        assert_eq!(gammafn(0.1), unsafe { c::gammafn(0.1) });
        assert_eq!(gammafn(1.0), unsafe { c::gammafn(1.0) });
    }

    #[test]
    fn test_lbeta() {
        assert_eq!(lbeta(0.0, 0.0), unsafe { c::lbeta(0.0, 0.0) });
        assert_eq!(lbeta(1.0, 1.0), unsafe { c::lbeta(1.0, 1.0) });
        assert_eq!(lbeta(10.0, 0.0), unsafe { c::lbeta(10.0, 0.0) });
        assert_eq!(lbeta(0.0, 10.0), unsafe { c::lbeta(0.0, 10.0) });
        assert_eq!(lbeta(1.0, 1.0), unsafe { c::lbeta(1.0, 1.0) });
        assert_eq!(lbeta(10.0, 10.0), unsafe { c::lbeta(10.0, 10.0) });
        assert_eq!(lbeta(100.0, 100.0), unsafe { c::lbeta(100.0, 100.0) });
        assert_eq!(lbeta(100.0, 1.0), unsafe { c::lbeta(100.0, 1.0) });
        assert_eq!(lbeta(1.0, 100.0), unsafe { c::lbeta(1.0, 100.0) });
        assert!(lbeta(-100.0, 0.0).is_nan());
        assert!(lbeta(0.0, -100.0).is_nan());
        assert!(lbeta(100.0, -100.0).is_nan());
    }

    #[test]
    fn test_lgammafn() {
        assert_eq!(lgammafn(0.0), unsafe { c::lgammafn(0.0) });
        assert_eq!(lgammafn(-1.0), unsafe { c::lgammafn(-1.0) });
        assert_eq!(lgammafn(1.0), unsafe { c::lgammafn(1.0) });
    }

    #[test]
    fn test_lgammafn_sign() {
        assert_eq!(lgammafn_sign(0.0, Some(&mut -1)), unsafe {
            c::lgammafn_sign(0.0, Some(&mut -1))
        });
        assert_eq!(lgammafn_sign(-1.0, Some(&mut -1)), unsafe {
            c::lgammafn_sign(-1.0, Some(&mut -1))
        });
        assert_eq!(lgammafn_sign(1.0, Some(&mut -1)), unsafe {
            c::lgammafn_sign(1.0, Some(&mut -1))
        });
        assert_eq!(lgammafn_sign(0.0, Some(&mut 1)), unsafe {
            c::lgammafn_sign(0.0, Some(&mut 1))
        });
        assert_eq!(lgammafn_sign(-1.0, Some(&mut 1)), unsafe {
            c::lgammafn_sign(-1.0, Some(&mut 1))
        });
        assert_eq!(lgammafn_sign(1.0, Some(&mut 1)), unsafe {
            c::lgammafn_sign(1.0, Some(&mut 1))
        });
    }

    #[test]
    fn test_log1pmx() {
        assert_eq!(log1pmx(0.0), unsafe { c::log1pmx(0.0) });
        assert_eq!(log1pmx(1.1), unsafe { c::log1pmx(1.1) });
        assert!(log1pmx(-1.23).is_nan());
        assert_eq!(log1pmx(1e-3), unsafe { c::log1pmx(1e-3) });
        // Test logcf via log1pmx.
        assert_eq!(log1pmx(0.91), unsafe { c::log1pmx(0.91) });
        assert_eq!(log1pmx(0.81), unsafe { c::log1pmx(0.81) });
    }

    #[test]
    fn test_pgamma() {
        assert!(pgamma(0.0, -1.0, 1.0, true, false).is_nan());
        assert!(pgamma(0.0, 1.0, -1.0, true, false).is_nan());
        assert!(pgamma(0.0, -1.0, -1.0, true, false).is_nan());
        assert_eq!(pgamma(0.1, 0.1, 1.0, false, false), unsafe {
            c::pgamma(0.1, 0.1, 1.0, 0, 0)
        });
        assert!(abs_diff_eq!(
            pgamma(0.65, 0.2, 0.34, false, false),
            unsafe { c::pgamma(0.65, 0.2, 0.34, 0, 0) },
            epsilon = 1e-15
        ));
        assert!(abs_diff_eq!(
            pgamma(3.21, 0.2, 0.34, false, false),
            unsafe { c::pgamma(3.21, 0.2, 0.34, 0, 0) },
            epsilon = 1e-15
        ));
        assert_eq!(pgamma(3.21, 0.2, 0.34, false, true), unsafe {
            c::pgamma(3.21, 0.2, 0.34, 0, 1)
        });
        assert_eq!(pgamma(123.0, 0.2, 0.34, false, true), unsafe {
            c::pgamma(123.0, 0.2, 0.34, 0, 1)
        });
    }

    #[test]
    fn test_pnorm() {
        assert_eq!(pnorm(0.0, 0.0, 1.0, true, false), 0.5);
        assert_eq!(pnorm(0.0, 0.0, 1.0, false, false), 0.5);
        assert_eq!(pnorm(0.0, 0.0, 1.0, false, false), unsafe {
            c::pnorm5(0.0, 0.0, 1.0, 0, 0)
        });
        assert_eq!(pnorm(0.65, 0.2, 0.34, false, false), unsafe {
            c::pnorm5(0.65, 0.2, 0.34, 0, 0)
        });
        assert_eq!(pnorm(3.21, 0.2, 0.34, false, false), unsafe {
            c::pnorm5(3.21, 0.2, 0.34, 0, 0)
        });
        assert_eq!(pnorm(3.21, 0.2, 0.34, false, true), unsafe {
            c::pnorm5(3.21, 0.2, 0.34, 0, 1)
        });
        assert_eq!(pnorm(123.0, 0.2, 0.34, false, true), unsafe {
            c::pnorm5(123.0, 0.2, 0.34, 0, 1)
        });
    }

    // TODO: add more tests
    #[test]
    fn test_pt() {
        // assert_eq!(pt(0.1, 1.0, false, false), unsafe { c::pt(0.1, 1.0, 0, 0) });
    }

    #[test]
    fn test_qnorm() {
        assert_eq!(qnorm(0.0, 0.5, 1.0, true, false), unsafe {
            c::qnorm5(0.0, 0.5, 1.0, 1, 0)
        });
        assert_eq!(qnorm(0.4, 0.5, 1.0, true, false), unsafe {
            c::qnorm5(0.4, 0.5, 1.0, 1, 0)
        });
        assert_eq!(qnorm(0.4, 0.5, 1.0, false, false), unsafe {
            c::qnorm5(0.4, 0.5, 1.0, 0, 0)
        });
        assert!(abs_diff_eq!(
            qnorm(-2.3, 0.5, 1.0, false, true),
            unsafe { c::qnorm5(-2.3, 0.5, 1.0, 0, 1) },
            epsilon = 1e-15
        ));
    }

    #[test]
    fn test_sinpi() {
        assert_eq!(sinpi(0.0), unsafe { c::sinpi(0.0) });
        assert_eq!(sinpi(0.234), unsafe { c::sinpi(0.234) });
    }

    #[test]
    fn test_tanpi() {
        assert_eq!(tanpi(0.0), unsafe { c::tanpi(0.0) });
        assert_eq!(tanpi(0.25), unsafe { c::tanpi(0.25) });
        assert_eq!(tanpi(0.234), unsafe { c::tanpi(0.234) });
    }
}
