use rmathlib::cospi;
use rmathlib::log1pmx;
use rmathlib::sinpi;
use rmathlib::tanpi;

mod c {
    extern "C" {
        pub fn log1pmx(x: f64) -> f64;
        pub fn cospi(x: f64) -> f64;
        pub fn sinpi(x: f64) -> f64;
        pub fn tanpi(x: f64) -> f64;
    }
}

#[test]
fn test_cospi() {
    assert_eq!(cospi(0.0), unsafe { c::cospi(0.0) });
    assert_eq!(cospi(0.234), unsafe { c::cospi(0.234) });
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
