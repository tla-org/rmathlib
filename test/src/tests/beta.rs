use rmathlib::lbeta;

mod c {
    extern "C" {
        pub fn lbeta(a: f64, b: f64) -> f64;
    }
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
