use rmathlib::bd0;
use rmathlib::chebyshev_eval;
use rmathlib::chebyshev_init;
use rmathlib::ebd0;
use rmathlib::rf_i1mach;
use rmathlib::stirlerr;

mod c {
    extern "C" {
        pub fn Rf_bd0(x: f64, np: f64) -> f64;
        pub fn Rf_chebyshev_eval(x: f64, a: *mut f64, n: i32) -> f64;
        pub fn Rf_chebyshev_init(dos: *mut f64, nos: i32, eta: f64) -> i32;
        pub fn Rf_ebd0(x: f64, M: f64, yh: *mut f64, yl: *mut f64);
        pub fn Rf_i1mach(i: i32) -> i32;
        pub fn Rf_stirlerr(n: f64) -> f64;
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
fn test_i1mach() {
    assert_eq!(rf_i1mach(1), unsafe { c::Rf_i1mach(1) });
    assert_eq!(rf_i1mach(2), unsafe { c::Rf_i1mach(2) });
    assert_eq!(rf_i1mach(3), unsafe { c::Rf_i1mach(3) });
    assert_eq!(rf_i1mach(4), unsafe { c::Rf_i1mach(4) });
    assert_eq!(rf_i1mach(5), unsafe { c::Rf_i1mach(5) });
    assert_eq!(rf_i1mach(6), unsafe { c::Rf_i1mach(6) });
    assert_eq!(rf_i1mach(7), unsafe { c::Rf_i1mach(7) });
    assert_eq!(rf_i1mach(8), unsafe { c::Rf_i1mach(8) });
    assert_eq!(rf_i1mach(9), unsafe { c::Rf_i1mach(9) });
    assert_eq!(rf_i1mach(10), unsafe { c::Rf_i1mach(10) });
    assert_eq!(rf_i1mach(11), unsafe { c::Rf_i1mach(11) });
    assert_eq!(rf_i1mach(12), unsafe { c::Rf_i1mach(12) });
    assert_eq!(rf_i1mach(13), unsafe { c::Rf_i1mach(13) });
    assert_eq!(rf_i1mach(14), unsafe { c::Rf_i1mach(14) });
    assert_eq!(rf_i1mach(15), unsafe { c::Rf_i1mach(15) });
    assert_eq!(rf_i1mach(16), unsafe { c::Rf_i1mach(16) });
}

#[test]
fn test_stirlerr() {
    assert_eq!(stirlerr(1.0), unsafe { c::Rf_stirlerr(1.0) });
    assert_eq!(stirlerr(2.0), unsafe { c::Rf_stirlerr(2.0) });
    assert_eq!(stirlerr(25.0), unsafe { c::Rf_stirlerr(25.0) });
    assert_eq!(stirlerr(50.0), unsafe { c::Rf_stirlerr(50.0) });
}
