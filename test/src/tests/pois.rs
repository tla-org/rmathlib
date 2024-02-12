use rmathlib::dpois;

mod c {
    extern "C" {
        pub fn dpois(x: f64, lambda: f64, give_log: bool) -> f64;
    }
}

#[test]
fn test_dpois() {
    assert!(dpois(-1.0, -1.0, false).is_nan());
    assert!(dpois(1.0, -1.0, false).is_nan());
    assert_eq!(dpois(1.0, 1.0, false), unsafe { c::dpois(1.0, 1.0, false) });
    assert_eq!(dpois(1.0, 1.0, true), unsafe { c::dpois(1.0, 1.0, true) });
}
