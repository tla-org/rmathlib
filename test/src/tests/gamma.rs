use approx::abs_diff_eq;
use rmathlib::gammafn;
use rmathlib::lgammacor;
use rmathlib::lgammafn;
use rmathlib::lgammafn_sign;
use rmathlib::pgamma;

mod c {
    extern "C" {
        pub fn Rf_lgammacor(x: f64) -> f64;
        pub fn gammafn(x: f64) -> f64;
        pub fn lgammafn(x: f64) -> f64;
        pub fn lgammafn_sign(x: f64, sgn: Option<&mut i32>) -> f64;
        pub fn pgamma(x: f64, alph: f64, scale: f64, lower_tail: i32, log_p: i32) -> f64;
    }
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
fn test_gammafn() {
    assert!(gammafn(-1.0).is_nan());
    assert!(gammafn(0.0).is_nan());
    assert_eq!(gammafn(0.1), unsafe { c::gammafn(0.1) });
    assert_eq!(gammafn(1.0), unsafe { c::gammafn(1.0) });
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
