#[cfg(test)]
mod test_math {

    mod c {
        extern "C" {
            pub fn cospi(x: f64) -> f64;
            pub fn sinpi(x: f64) -> f64;
            pub fn tanpi(x: f64) -> f64;
        }
    }

    #[test]
    fn test_cospi() {
        assert_eq!(rmathlib::cospi(0.0), unsafe { c::cospi(0.0) });
        assert_eq!(rmathlib::cospi(0.234), unsafe { c::cospi(0.234) });
        assert_eq!(rmathlib::sinpi(0.0), unsafe { c::sinpi(0.0) });
        assert_eq!(rmathlib::sinpi(0.234), unsafe { c::sinpi(0.234) });
        assert_eq!(rmathlib::tanpi(0.0), unsafe { c::tanpi(0.0) });
        assert_eq!(rmathlib::tanpi(0.25), unsafe { c::tanpi(0.25) });
        assert_eq!(rmathlib::tanpi(0.234), unsafe { c::tanpi(0.234) });
    }
}
