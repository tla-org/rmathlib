use approx::abs_diff_eq;
use rmathlib::dnorm;
use rmathlib::pnorm;
use rmathlib::qnorm;

mod c {
    extern "C" {
        pub fn dnorm4(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64;
        pub fn pnorm5(x: f64, mu: f64, sigma: f64, lower_tail: i32, log_p: i32) -> f64;
        pub fn qnorm5(p: f64, mu: f64, sigma: f64, lower_tail: i32, log_p: i32) -> f64;
    }
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
