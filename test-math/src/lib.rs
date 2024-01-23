#[cfg(test)]
mod test_math {

    extern "C" {
        fn cospi(x: f64) -> f64;
        fn sinpi(x: f64) -> f64;
        fn tanpi(x: f64) -> f64;
    }

    #[test]
    fn test_cospi() {
        assert_eq!(math::cospi(0.0), unsafe { cospi(0.0) });
        assert_eq!(math::cospi(0.234), unsafe { cospi(0.234) });
        assert_eq!(math::sinpi(0.0), unsafe { sinpi(0.0) });
        assert_eq!(math::sinpi(0.234), unsafe { sinpi(0.234) });
        assert_eq!(math::tanpi(0.0), unsafe { tanpi(0.0) });
        assert_eq!(math::tanpi(0.25), unsafe { tanpi(0.25) });
        assert_eq!(math::tanpi(0.234), unsafe { tanpi(0.234) });
    }
}
