#[cfg(test)]
mod test_math {
    mod c {
        extern "C" {
            pub fn pbeta(x: f64, a: f64, b: f64, lower_tail: i32, log_p: i32) -> f64;
        }
    }

    use approx::assert_abs_diff_eq;
    use rmathlib::*;
    use statrs::statistics::Statistics;

    #[test]
    // This test is based on `toms708_test.f` subroutine `TEST01`.
    // Test BRATIO computes the Beta ratio function.
    fn subroutine_test01() {
        let a: f64 = 5.3;
        let b: f64 = 10.1;

        for i in 1..=50 {
            let x: f64 = (i as f64) / 100.0;
            let y: f64 = 0.5 + (0.5 - x);
            let log_p = false;

            let mut _w: f64 = 0.0;
            let mut _w1: f64 = 0.0;
            let mut ierr: i32 = 1;
            bratio(a, b, x, y, &mut _w, &mut _w1, &mut ierr, log_p);
            assert!(ierr == 0);
        }
    }

    #[test]
    fn test_pbeta_and_toms708() {
        // Manual values obtained from R 4.3.2.
        assert_eq!(pbeta(0.9, 0.1, 0.1, false, false), 0.40638509393627587);
        // assert!(pbeta(0.9, 0.0, 1.0, false, false).is_nan());
        // assert_eq!(pbeta(0.1, 0.5, 0.5, true, true), 0.1285223);

        fn helper(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) {
            let rs = pbeta(x, a, b, lower_tail, log_p);
            let c = unsafe { c::pbeta(x, a, b, lower_tail as i32, log_p as i32) };
            // If both NaN, then avoid comparing since that would fail.
            if !(rs.is_nan() && c.is_nan()) {
                assert_eq!(
                    rs, c,
                    "pbeta({}, {}, {}, {}, {})",
                    x, a, b, lower_tail, log_p
                );
            }
        }

        // These functions do not call `bratio`.
        helper(0.9, 0.0, 1.0, false, false);
        helper(1.1, 0.0, 1.0, false, false);
        helper(0.0, -1.0, 1.0, false, false);
        helper(f64::NAN, 1.0, 1.0, false, false);
        helper(0.1, 0.0, 1.0, false, false);

        // These functions call `bratio`.
        // Thanks to the `pbeta` definition in `pbeta.rs`, we can compare the
        // outcome of `bratio` against the outcome of `pbeta` in R.
        // The values are obtained via R version 4.3.4.
        //
        // `pbeta.rs` is mostly handling edge cases and then calling `bratio`.
        // Cases which end up in `bratio` are where 0 < a < Inf and 0 < b < Inf.

        // R> pbeta(0.5, 1.0, 1.0, lower.tail = TRUE, log.p = FALSE)
        assert_eq!(pbeta(0.5, 1.0, 1.0, true, false), 0.5);

        let epsilon = 1e-12;
        assert_abs_diff_eq!(
            pbeta(0.01, 0.01, 0.01, true, false),
            // R> sprintf("%.13f", pbeta(0.01, 0.01, 0.01, lower.tail = TRUE, log.p = FALSE))
            0.4776207614162,
            epsilon = epsilon
        );

        assert_abs_diff_eq!(
            pbeta(0.001, 0.001, 0.01, true, false),
            // R> sprintf("%.13f", pbeta(0.001, 0.001, 0.01, lower.tail = TRUE, log.p = FALSE))
            0.9028483975306,
            epsilon = epsilon
        );

        assert_abs_diff_eq!(
            pbeta(0.1, 0.8, 2.0, true, true),
            // R> sprintf("%.13f", pbeta(0.1, 0.8, 2.0, lower.tail=TRUE, log.p=TRUE))
            -1.2997437835699,
            epsilon = epsilon
        );

        assert_abs_diff_eq!(
            pbeta(0.1, 0.8, 2.0, false, true),
            // R> sprintf("%.13f", pbeta(0.1, 0.8, 2.0, lower.tail=FALSE, log.p=TRUE))
            -0.3182809860569,
            epsilon = epsilon
        );

        // Based on a test in `d-p-q-r-tst-2.R` from the R source code.
        for x in [0.01, 0.10, 0.25, 0.40, 0.55, 0.71, 0.98] {
            assert_abs_diff_eq!(
                pbeta(x, 0.8, 2.0, false, true),
                pbeta(1.0 - x, 2.0, 0.8, true, true),
                epsilon = epsilon,
            )
        }

        // Based on a test in `d-p-q-r-tst-2.R` from the R source code.
        assert_abs_diff_eq!(
            pbeta(256.0 / 1024.0, 3.0, 2200.0, false, true),
            // R> sprintf("%.13f", pbeta(256/1024, 3, 2200, lower.tail=FALSE, log.p=TRUE))
            -620.9697808693397,
            epsilon = epsilon
        );

        assert_abs_diff_eq!(
            pbeta(512.0 / 1024.0, 3.0, 2200.0, false, true),
            // R> sprintf("%.13f", pbeta(512/1024, 3, 2200, lower.tail=FALSE, log.p=TRUE))
            -1511.6085416971891,
            epsilon = epsilon
        );

        assert_abs_diff_eq!(
            pbeta(768.0 / 1024.0, 3.0, 2200.0, false, true),
            // R> sprintf("%.13f", pbeta(768/1024, 3, 2200, lower.tail=FALSE, log.p=TRUE))
            -3035.7220144978146,
            epsilon = epsilon
        );

        assert!(pbeta(1024.0 / 1024.0, 3.0, 2200.0, false, true).is_infinite());

        fn diff(x: Vec<f64>) -> Vec<f64> {
            let mut result = Vec::new();
            let n = x.len();
            for i in 0..n - 1 {
                result.push(x[i + 1] - x[i]);
            }
            result
        }

        fn log(x: Vec<f64>) -> Vec<f64> {
            x.iter().map(|x| x.ln()).collect()
        }

        // Based on a test in `d-p-q-r-tst-2.R` at line 330 from the R source code.
        // pbeta(x, a, b, log=TRUE) for small x and a is ~ log-linear.
        let x = (10..=200).map(|n| 2.0_f64.powf(-n as f64));
        for a in [1e-8, 1e-12, 16e-16, 4e-16] {
            for b in [0.6, 1.0, 2.0, 10.0] {
                let xs = x.clone().map(|x| pbeta(x, a, b, true, true));
                let dp = diff(xs.collect());
                let sd = dp.clone().population_std_dev();
                let m = dp.mean();
                assert!(sd / m < 0.0007);
            }
        }

        // Based on a test in `d-p-q-r-tst-2.R` at line 401 from the R source code.
        // pbeta(x, <small a>,<small b>, .., log):
        let a = (1..=25).map(|i| 2_f64.powi(-(90 + i)));
        let b = 2_f64.powi(-60);
        let pb = a.clone().map(|a| pbeta(0.5, a, b, true, true));
        let ldp = diff(log(diff(pb.collect())));
        for l in ldp {
            let val = l - (1_f64 / 2_f64).ln();
            assert!(val.abs() < 1e-9)
        }

        // Based on a test in `d-p-q-r-tst-2.R` at line 455 from the R source code.
        // qbeta(*, a,b) when  a,b << 1 : can easily fail
        {
            let q1 = 9.094955e-13;
            assert_abs_diff_eq!(
                pbeta(q1, 0.125, 2f64.powi(-26), true, false),
                2f64.powi(-28),
                epsilon = epsilon
            )
        }

        // Some additional tests.
        helper(0.1, 0.5, 0.5, false, false);
        helper(0.1, 0.5, 0.5, true, false);
        helper(0.1, 0.5, 0.5, true, true);
        helper(0.1, 1.001, 0.99, true, true);
        helper(0.5, 10000.0, 0.2, true, true);
        helper(0.5, 10000.0, 0.2, true, false);
        helper(0.5, 10000.0, 0.2, false, false);
        helper(0.0, 0.0, 0.2, false, false);

         // Based on a test in `d-p-q-r-tst-2.R` at line 850 from the R source code.
        assert_eq!(pbeta(0.0, 0.0, 3.0, true, false), 1.0);
        assert_eq!(pbeta(1.0, 0.1, 0.0, true, false), 1.0);
        assert_eq!(pbeta(1.1, 3.0, 0.0, true, false), 1.0);
        assert_eq!(pbeta(0.0, 0.0, 0.0, true, false), 0.5);
        assert_eq!(pbeta(1.0, 0.0, 0.0, true, false), 1.0);
    }
}
