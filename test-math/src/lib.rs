#[cfg(test)]
mod test_math {

    extern "C" {
        fn cospi(x: f64) -> f64;
    }

    #[test]
    fn test_cospi() {
        assert_eq!(math::cospi(0.0), unsafe { cospi(0.0) });
        assert_eq!(math::cospi(0.234), unsafe { cospi(0.234) });
    }
}
