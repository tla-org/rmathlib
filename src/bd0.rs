//!
//! AUTHORS
//! Catherine Loader, catherine@research.bell-labs.com, October 23, 2000. [ bd0() ]
//! Morten Welinder, see Bugzilla PR#15628, 2014                          [ebd0() ]
//!
//! Merge in to R (and much more):
//!
//! Copyright (C) 2000-2022 The R Core Team
//!
//! This program is free software; you can redistribute it and/or modify
//! it under the terms of the GNU General Public License as published by
//! the Free Software Foundation; either version 2 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU General Public License for more details.
//!
//! You should have received a copy of the GNU General Public License
//! along with this program; if not, a copy is available at
//! https://www.R-project.org/Licenses/

use crate::libc::*;
use crate::nmath::*;
use crate::pgamma::log1pmx;
use crate::rmath::*;
use libm::frexp;
use libm::ldexp;

/// Calculates a stable deviance part using a method that reduces relative error
///
/// Evaluates the "deviance part"
/// bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
///           =  x * log(x/M) + M - x
/// where M = E[X] = n*p (or = lambda), for x, M > 0
///
/// in a manner that should be stable (with small relative error)
/// for all x and M=np. In particular for x/np close to 1, direct
/// evaluation fails, and evaluation is based on the Taylor series
/// of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
///
/// Martyn Plummer had the nice idea to use log1p() and Martin Maechler
/// emphasized the extra need to control cancellation.
///
/// MP:   t := (x-M)/M  ( <==> 1+t = x/M  ==>
///
/// bd0 = M*[ x/M * log(x/M) + 1 - (x/M) ] = M*[ (1+t)*log1p(t) + 1 - (1+t) ]
///     = M*[ (1+t)*log1p(t) - t ] =: M * p1log1pm(t) =: M * p1l1(t)
/// MM: The above is very nice, as the "simple" p1l1() function would be useful
///   to have available in a fast numerical stable way more generally.
pub fn bd0(x: f64, np: f64) -> f64 {
    if !r_finite(x) || !r_finite(np) || np == 0.0 {
        return ml_warn_return_nan();
    }

    if fabs(x - np) < 0.1 * (x + np) {
        let mut v: f64 = (x - np) / (x + np);
        let mut s: f64 = (x - np) * v;
        if fabs(s) < DBL_MIN {
            return s;
        }
        let mut ej: f64 = 2.0 * x * v;
        v *= v; // v^2
        for j in 1..1000 {
            // Taylor series; 1000: no infinite loop
            // as |v| < 0.1, v^2000 is "zero".
            ej *= v;
            let s_: f64 = s;
            s += ej / ((j << 1) + 1) as f64;
            if s == s_ {
                // Last term was effectively 0.
                return s;
            }
        }
    }
    println!("bd0: T.series failed to converge in 1000 iterations");
    x * log(x / np) + np - x
}

#[allow(clippy::approx_constant)]
const BD0_SCALE: [[f32; 4]; 129] = [
    [
        0.693_147_2,
        -1.904_654_2e-9,
        -8.783_184e-17,
        3.061_840_7e-24,
    ], /* 128: log(2048/1024.) */
    [
        0.685_304_05,
        -4.257_826_4e-8,
        -1.172_310_5e-15,
        6.203_392_6e-23,
    ], /* 129: log(2032/1024.) */
    [0.677_398_8, 2.274_189e-8, 1.441_192e-15, 7.046_384_5e-23], /* 130: log(2016/1024.) */
    [
        0.669_930_6,
        -4.829_385_6e-8,
        -8.664_795_5e-16,
        7.049_558e-24,
    ], /* 131: log(2001/1024.) */
    [
        0.662_406_1,
        -4.791_602_5e-8,
        -2.161_508_2e-15,
        7.092_968_4e-23,
    ], /* 132: log(1986/1024.) */
    [0.654_824_5, 6.237_715_2e-9, 1.180_67e-16, 6.660_323_5e-25], /* 133: log(1971/1024.) */
    [
        0.647_185_1,
        -4.220_867e-8,
        -1.381_758_9e-15,
        -1.115_993_24e-23,
    ], /* 134: log(1956/1024.) */
    [0.640_001_9, -4.979_17e-8, 4.887_014e-16, 1.684_746_6e-23], /* 135: log(1942/1024.) */
    [
        0.632_766_7,
        -5.406_177_2e-8,
        -2.722_454_5e-15,
        -2.907_078e-23,
    ], /* 136: log(1928/1024.) */
    [0.624_956_13, 3.285_935_4e-8, 1.201_673e-16, -8.668_233e-24], /* 137: log(1913/1024.) */
    [
        0.618_137_36,
        -6.406_189_5e-11,
        3.550_581e-18,
        -3.623_679_4e-25,
    ], /* 138: log(1900/1024.) */
    [0.610_741_6, 4.229_853_8e-8, 1.548_143_3e-15, -5.446_353e-23], /* 139: log(1886/1024.) */
    [
        0.603_290_8,
        5.515_817_5e-8,
        2.119_363_7e-15,
        1.247_109_8e-22,
    ], /* 140: log(1872/1024.) */
    [0.596_322_2, 3.281_353_9e-9, -1.403_347e-16, 5.544_062e-24], /* 141: log(1859/1024.) */
    [
        0.589_304_57,
        4.307_998_6e-8,
        -3.244_665_2e-15,
        -2.569_782_4e-23,
    ], /* 142: log(1846/1024.) */
    [0.582_237_5, -3.983_067_4e-8, 2.938_790_6e-15, 1.612_443e-22], /* 143: log(1833/1024.) */
    [0.575_12, 2.242_384e-9, -2.204_516_4e-17, 4.078_543_8e-25], /* 144: log(1820/1024.) */
    [0.568_504_7, 4.422_870_6e-8, 2.787_995_7e-16, -9.741_709e-24], /* 145: log(1808/1024.) */
    [
        0.561_845_4,
        2.147_161_4e-8,
        1.337_491_9e-15,
        -2.326_922_3e-23,
    ], /* 146: log(1796/1024.) */
    [0.554_580_8, 4.577_835e-9, 2.633_146e-16, 2.070_995_9e-23], /* 147: log(1783/1024.) */
    [0.547_827_84, -7.668e-9, 6.199_095_3e-16, 6.341_979e-24],   /* 148: log(1771/1024.) */
    [0.541_029, -3.730_232e-8, -3.380_178_1e-15, 1.446_919_8e-22], /* 149: log(1759/1024.) */
    [0.534_755_7, 4.382_892e-8, -8.601_821e-16, -8.775_950_6e-24], /* 150: log(1748/1024.) */
    [0.527_867_1, 1.083_971_5e-8, -4.480_281e-16, 3.884_099_6e-23], /* 151: log(1736/1024.) */
    [0.521_510_5, 4.203_159_4e-8, 1.211_01e-15, -2.320_674_5e-23], /* 152: log(1725/1024.) */
    [
        0.514_529_7,
        -1.232_294_1e-8,
        2.820_084_5e-16,
        1.880_026_6e-23,
    ], /* 153: log(1713/1024.) */
    [
        0.508_087_5,
        -1.229_712_7e-8,
        -4.295_835_6e-16,
        -1.803_662_6e-23,
    ], /* 154: log(1702/1024.) */
    [0.501_603_5, 5.914_538e-8, 1.272_856_1e-15, -5.824_094_4e-23], /* 155: log(1691/1024.) */
    [
        0.495_077_25,
        1.440_985_1e-8,
        -6.381_91e-17,
        -3.500_509_4e-25,
    ], /* 156: log(1680/1024.) */
    [
        0.489_107_07,
        2.745_798_6e-8,
        -1.447_041_9e-15,
        4.211_867e-23,
    ], /* 157: log(1670/1024.) */
    [
        0.482_498_47,
        1.762_245_5e-8,
        4.128_675_2e-16,
        2.608_269_2e-23,
    ], /* 158: log(1659/1024.) */
    [
        0.476_452_53,
        -1.247_024_35e-8,
        -9.162_194e-17,
        3.765_782_5e-24,
    ], /* 159: log(1649/1024.) */
    [0.469_759_46, -5.450_354e-9, -4.304_847e-16, 2.074_710_3e-25], /* 160: log(1638/1024.) */
    [
        0.463_635_74,
        -1.701_304_7e-9,
        7.601_623e-18,
        -2.041_587_9e-25,
    ], /* 161: log(1628/1024.) */
    [
        0.457_474_3,
        6.943_684_5e-10,
        -2.546_131e-17,
        5.541_253_3e-25,
    ], /* 162: log(1618/1024.) */
    [0.451_274_63, 1.073_186_5e-8, 6.374e-16, 3.901_547_3e-23],  /* 163: log(1608/1024.) */
    [
        0.445_036_3,
        2.865_065_7e-8,
        -9.155_353e-16,
        -4.736_587_8e-23,
    ], /* 164: log(1598/1024.) */
    [
        0.439_388_33,
        2.618_613_2e-8,
        1.361_960_2e-15,
        -5.267_279_5e-23,
    ], /* 165: log(1589/1024.) */
    [
        0.433_075_2,
        1.906_573_4e-8,
        1.014_319_2e-15,
        1.014_567_1e-22,
    ], /* 166: log(1579/1024.) */
    [
        0.427_359_1,
        -1.141_359_4e-8,
        1.324_204_5e-16,
        -1.024_005_6e-23,
    ], /* 167: log(1570/1024.) */
    [
        0.420_969_3,
        -1.277_850_9e-8,
        6.143_525_7e-16,
        1.419_242_2e-23,
    ], /* 168: log(1560/1024.) */
    [0.415_183_37, -7.767_916e-9, 5.955_443e-16, 2.732_668e-23], /* 169: log(1551/1024.) */
    [
        0.409_363_75,
        1.880_755e-9,
        1.915_333_1e-16,
        -5.620_806_3e-24,
    ], /* 170: log(1542/1024.) */
    [0.403_510_1, -2.041_660_4e-8, -2.934_05e-16, 1.894_347e-24], /* 171: log(1533/1024.) */
    [0.397_621_93, 1.001_600_1e-9, 2.286_324_2e-17, 9.458_134e-25], /* 172: log(1524/1024.) */
    [
        0.391_698_9,
        1.545_909_7e-8,
        1.096_282_3e-15,
        3.108_302_2e-23,
    ], /* 173: log(1515/1024.) */
    [
        0.386_404_4,
        -2.076_412_4e-9,
        1.507_346_5e-16,
        7.412_449_5e-24,
    ], /* 174: log(1507/1024.) */
    [
        0.380_414_37,
        -8.244_395e-9,
        1.486_622_5e-16,
        -3.927_292_7e-24,
    ], /* 175: log(1498/1024.) */
    [0.374_388_22, 9.158_53e-9, 5.656_919e-16, 3.421_347_5e-23], /* 176: log(1489/1024.) */
    [
        0.369_001_03,
        -2.225_359_1e-8,
        6.231_405_4e-16,
        -2.156_475_1e-23,
    ], /* 177: log(1481/1024.) */
    [
        0.363_584_64,
        -2.677_872_9e-8,
        -9.943_908e-16,
        -4.704_929_7e-24,
    ], /* 178: log(1473/1024.) */
    [
        0.357_455_9,
        -2.033_036_1e-8,
        -1.579_449_2e-15,
        6.318_678e-23,
    ], /* 179: log(1464/1024.) */
    [0.351_976_4, 2.850_385_9e-8, -9.566_435e-16, -6.409_595e-24], /* 180: log(1456/1024.) */
    [
        0.346_466_78,
        -1.236_265_3e-8,
        -6.003_368e-16,
        -2.860_901_5e-24,
    ], /* 181: log(1448/1024.) */
    [
        0.340_926_6,
        -6.110_413_3e-10,
        1.746_713_6e-17,
        1.996_258_7e-25,
    ], /* 182: log(1440/1024.) */
    [
        0.335_355_52,
        2.167_272_6e-8,
        -1.091_877_3e-15,
        -3.047_574_8e-23,
    ], /* 183: log(1432/1024.) */
    [0.330_455_3, -1.608_884e-8, -3.833_435e-16, -7.683_742e-24], /* 184: log(1425/1024.) */
    [
        0.324_825_4,
        2.801_670_4e-8,
        -2.072_572_1e-16,
        1.316_077_8e-23,
    ], /* 185: log(1417/1024.) */
    [0.319_163_68, 2.622_263e-8, -1.399_522_3e-15, 8.599_884e-23], /* 186: log(1409/1024.) */
    [
        0.314_183_24,
        2.682_662_6e-8,
        -9.792_556e-16,
        2.295_496_1e-23,
    ], /* 187: log(1402/1024.) */
    [
        0.308_460_77,
        1.368_350_9e-8,
        5.591_995e-16,
        -1.193_870_14e-23,
    ], /* 188: log(1394/1024.) */
    [0.303_426_62, -8.629_042e-9, -5.222_554e-16, 3.228_770_8e-23], /* 189: log(1387/1024.) */
    [
        0.298_366_96,
        8.688_424e-9,
        2.764_116_8e-16,
        -1.017_185_9e-23,
    ], /* 190: log(1380/1024.) */
    [
        0.292_553,
        -4.916_314_5e-9,
        2.562_284_7e-16,
        -2.634_157_6e-23,
    ], /* 191: log(1372/1024.) */
    [
        0.287_437_92,
        -1.378_239_6e-8,
        7.290_393_5e-16,
        -4.431_978e-24,
    ], /* 192: log(1365/1024.) */
    [0.282_296_48, 2.377_086_6e-8, 6.449_228_4e-16, 3.417_537e-23], /* 193: log(1358/1024.) */
    [0.277_128_52, 1.473_302_9e-8, 6.793_364e-16, 3.593_898e-24], /* 194: log(1351/1024.) */
    [0.271_933_73, -1.893_332e-8, 7.779_394e-16, 2.591_972_4e-23], /* 195: log(1344/1024.) */
    [
        0.266_711_77,
        1.430_038_7e-11,
        -7.458_877e-19,
        4.247_418_8e-26,
    ], /* 196: log(1337/1024.) */
    [0.262_214, 7.802_219e-9, 5.022_431e-16, -2.406_317_5e-23],  /* 197: log(1331/1024.) */
    [0.256_940_9, 2.961_805e-8, 7.279_528e-16, 2.638_556_1e-23], /* 198: log(1324/1024.) */
    [0.251_639_9, -6.447_878e-9, 4.122_028_1e-16, 7.427_559e-24], /* 199: log(1317/1024.) */
    [0.247_073_68, -1.998_183e-9, -1.090_975_8e-16, 6.236_552e-24], /* 200: log(1311/1024.) */
    [0.241_719_93, 5.523_086e-9, -2.686_547_6e-16, -2.764_496e-24], /* 201: log(1304/1024.) */
    [
        0.237_108_08,
        1.008_537_4e-8,
        -4.775_627e-16,
        2.330_620_6e-23,
    ], /* 202: log(1298/1024.) */
    [
        0.231_700_6,
        -1.394_638_4e-8,
        1.097_092_1e-17,
        3.229_909_4e-25,
    ], /* 203: log(1291/1024.) */
    [
        0.227_042_2,
        -6.451_285e-9,
        -4.252_994_8e-16,
        -1.026_025_5e-23,
    ], /* 204: log(1285/1024.) */
    [0.222_361_98, 1.411_064_5e-8, 6.025_569e-16, 1.938_518e-23], /* 205: log(1279/1024.) */
    [0.217_659_8, -8.286_783e-9, 5.232_364e-16, 5.078_443_5e-23], /* 206: log(1273/1024.) */
    [0.212_145_8, -8.254_219e-9, 3.255_553_1e-16, 1.571_43e-23], /* 207: log(1266/1024.) */
    [0.207_395_2, -1.614_928e-9, 2.113_159_3e-17, 1.427_561_7e-24], /* 208: log(1260/1024.) */
    [0.202_621_9, 8.597_64e-9, -3.380_462e-16, 2.562_323_6e-24], /* 209: log(1254/1024.) */
    [
        0.197_825_73,
        1.348_296_6e-8,
        -3.202_457e-16,
        -2.571_225_2e-23,
    ], /* 210: log(1248/1024.) */
    [
        0.193_006_46,
        9.956_86e-10,
        9.001_674_6e-17,
        -3.754_797_7e-24,
    ], /* 211: log(1242/1024.) */
    [
        0.188_972_56,
        4.241_536e-9,
        3.808_681_5e-16,
        -2.114_740_3e-23,
    ], /* 212: log(1237/1024.) */
    [
        0.184_110_31,
        6.931_055e-9,
        -3.478_485_9e-16,
        2.466_594_3e-23,
    ], /* 213: log(1231/1024.) */
    [
        0.179_224_31,
        5.073_923_5e-9,
        3.222_133e-16,
        -1.037_900_9e-23,
    ], /* 214: log(1225/1024.) */
    [0.174_314_32, 3.794_385_3e-9, 3.190_067e-16, 2.029_271_5e-23], /* 215: log(1219/1024.) */
    [
        0.170_204_16,
        3.422_334_4e-9,
        -1.884_641_7e-16,
        1.141_531_5e-23,
    ], /* 216: log(1214/1024.) */
    [
        0.165_249_59,
        -1.321_003_9e-8,
        -2.321_395_4e-16,
        3.043_054_2e-24,
    ], /* 217: log(1208/1024.) */
    [0.160_270_3, 6.007_922e-9, -7.521_048e-17, -1.264_910_6e-25], /* 218: log(1202/1024.) */
    [
        0.156_101_91,
        -1.230_153_6e-8,
        3.017_561_8e-16,
        -8.633_806_5e-24,
    ], /* 219: log(1197/1024.) */
    [
        0.151_916_06,
        -1.484_557_2e-8,
        -3.265_830_3e-16,
        -1.526_815_2e-23,
    ], /* 220: log(1192/1024.) */
    [
        0.146_869_78,
        -4.674_899_6e-9,
        -4.294_291_3e-16,
        1.328_296e-23,
    ], /* 221: log(1186/1024.) */
    [0.142_645, 9.186_472e-9, -8.049_373e-16, 1.437_998_8e-23],  /* 222: log(1181/1024.) */
    [0.138_402_31, 9.865_117e-9, -8.837_306_5e-16, 7.295_325e-24], /* 223: log(1176/1024.) */
    [0.133_287_22, 9.990_351e-10, 3.316_946_6e-17, 2.735_144e-24], /* 224: log(1170/1024.) */
    [
        0.129_004_57,
        -7.461_208_5e-9,
        -6.212_113e-16,
        1.855_187_3e-24,
    ], /* 225: log(1165/1024.) */
    [
        0.124_703_48,
        -3.292_446_3e-9,
        -7.404_12e-17,
        1.324_695_6e-24,
    ], /* 226: log(1160/1024.) */
    [
        0.120_383_814,
        3.379_199_1e-9,
        1.621_498_2e-16,
        -6.000_706_7e-24,
    ], /* 227: log(1155/1024.) */
    [
        0.116_045_415,
        3.563_839_2e-10,
        -7.354_22e-18,
        7.794_312_4e-26,
    ], /* 228: log(1150/1024.) */
    [0.111_688_11, 3.136_766e-9, -8.994_406e-19, -7.920_937_6e-26], /* 229: log(1145/1024.) */
    [0.107_311_74, -4.728_528e-9, -4.297_634e-16, 5.351_143_3e-24], /* 230: log(1140/1024.) */
    [
        0.102_916_12,
        2.833_200_8e-9,
        4.925_742_8e-17,
        2.794_436_8e-24,
    ], /* 231: log(1135/1024.) */
    [0.098_501_1, 4.970_725e-9, 4.130_512_8e-16, 6.313_449_6e-25], /* 232: log(1130/1024.) */
    [
        0.094_066_516,
        -6.525_851e-9,
        -1.349_281_7e-16,
        -9.079_65e-24,
    ], /* 233: log(1125/1024.) */
    [
        0.089_612_156,
        2.536_961_8e-9,
        1.611_066_5e-16,
        -5.189_504_5e-24,
    ], /* 234: log(1120/1024.) */
    [
        0.086_034_34,
        -5.304_795_7e-9,
        5.127_575_5e-17,
        1.463_615_5e-24,
    ], /* 235: log(1116/1024.) */
    [
        0.081_543_98,
        2.011_215_6e-9,
        8.065_769_3e-17,
        -3.015_032e-24,
    ], /* 236: log(1111/1024.) */
    [0.077_033_37, 5.749_566e-9, -2.503_851e-16, -1.846_143_1e-23], /* 237: log(1106/1024.) */
    [
        0.072_502_33,
        1.177_662_6e-9,
        -3.525_477e-17,
        1.316_407_8e-24,
    ], /* 238: log(1101/1024.) */
    [
        0.068_862_66,
        -7.043_545_3e-9,
        2.497_124e-16,
        1.068_688_25e-23,
    ], /* 239: log(1097/1024.) */
    [0.064_294_35, -2.422_082e-9, -2.055_589_6e-16, 8.602_908e-24], /* 240: log(1092/1024.) */
    [
        0.060_624_62,
        7.905_943e-12,
        -8.270_443_3e-19,
        -2.353_382_1e-26,
    ], /* 241: log(1088/1024.) */
    [
        0.056_018_442,
        -5.139_745_3e-10,
        -3.811_651_6e-17,
        1.907_219_5e-24,
    ], /* 242: log(1083/1024.) */
    [
        0.052_318_163,
        -3.557_432_6e-9,
        9.191_156e-17,
        -5.321_464e-24,
    ], /* 243: log(1079/1024.) */
    [
        0.047_673_47,
        -1.802_634_9e-9,
        1.032_963_4e-16,
        -2.228_357e-24,
    ], /* 244: log(1074/1024.) */
    [
        0.043_942_124,
        -1.795_005_7e-9,
        -5.381_740_3e-17,
        -1.399_619_7e-24,
    ], /* 245: log(1070/1024.) */
    [0.040_196_8, 3.845_193e-10, -2.448_545_3e-17, -7.386_769e-26], /* 246: log(1066/1024.) */
    [
        0.035_495_333,
        -3.590_165_4e-10,
        -2.073_207_8e-17,
        -2.412_097_2e-26,
    ], /* 247: log(1061/1024.) */
    [
        0.031_718_18,
        6.872_350_5e-10,
        -6.430_478e-18,
        1.350_869_2e-25,
    ], /* 248: log(1057/1024.) */
    [
        0.027_926_706,
        7.568_773_4e-10,
        -4.220_316e-17,
        2.534_760_6e-24,
    ], /* 249: log(1053/1024.) */
    [
        0.024_120_804,
        -1.125_570_7e-9,
        4.897_006e-17,
        1.417_221_5e-24,
    ], /* 250: log(1049/1024.) */
    [
        0.019_342_963,
        1.906_861e-10,
        -1.063_594_6e-17,
        -5.300_489_5e-25,
    ], /* 251: log(1044/1024.) */
    [
        0.015_504_187,
        -4.370_104_8e-10,
        6.611_061_5e-18,
        2.539_808_7e-25,
    ], /* 252: log(1040/1024.) */
    [
        0.011_650_616,
        9.168_89e-10,
        -1.584_869_8e-17,
        -1.350_491_6e-24,
    ], /* 253: log(1036/1024.) */
    [
        0.007_782_140_7,
        -3.046_577_4e-10,
        7.793_436e-18,
        4.660_1e-25,
    ], /* 254: log(1032/1024.) */
    [
        0.003_898_640_6,
        -2.132_467_8e-10,
        1.254_165_8e-19,
        8.745_035_4e-27,
    ], /* 255: log(1028/1024.) */
    [0.0, 0.0, 0.0, 0.0],                                        /* log(1024/1024) = log(1) = 0 */
];

fn add1(d_: f64, yh: &mut f64, yl: &mut f64) {
    let d = d_;
    let d1 = (d + 0.5).floor();
    let d2 = d - d1; // In [-0.5, 0.5)

    *yh += d1;
    *yl += d2;
}

/// Compute x * log (x / M) + (M - x) aka -x * log1pmx ((M - x) / x)
///
/// Returns `(yh, yl)`.
///
/// Unlike the C version, this function does not take mutable `yl` and `yh` arguments.
/// This is to make the interface easier to use.
pub fn ebd0(x: f64, m: f64) -> (f64, f64) {
    let sb: i32 = 10;
    let s: f64 = (1u32 << sb) as f64;
    let n: i32 = 128;

    let mut yh = 0.0;
    let mut yl = 0.0;

    if x == m {
        return (yh, yl);
    }

    if x == 0.0 {
        yh = m;
        return (yh, yl);
    }

    if m == 0.0 {
        yh = ML_POSINF;
        return (yh, yl);
    }

    // NB: m/x overflow handled above; underflow should be handled by fg = Inf.
    let (r, e) = frexp(m / x);

    // Prevent later overflow.
    if (M_LN2 * (-e as f64)) > 1.0 + DBL_MAX / x {
        yh = ML_POSINF;
        return (yh, yl);
    }

    let i: i32 = ((r - 0.5) * (2 * n) as f64 + 0.5).floor() as i32;
    // Now, 0 <= i <= n.
    let f: f64 = (s / (0.5 + i as f64 / (2.0 * n as f64) + 0.5)).floor();
    let fg: f64 = ldexp(f, -(e + sb)); // ldexp(f, E) := f * 2^E.
    if fg == ML_POSINF {
        yh = fg;
        return (yh, yl);
    }

    // We now have (M * fg / x) close to 1.  */
    //
    // We need to compute this:
    // (x/M)^x * exp(M-x) =
    // (M/x)^-x * exp(M-x) =
    // (M*fg/x)^-x * (fg)^x * exp(M-x) =
    // (M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)
    //
    // In log terms:
    // log((x/M)^x * exp(M-x)) =
    // log((M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)) =
    // log((M*fg/x)^-x * exp(M*fg-x)) + x*log(fg) + (M-M*fg) =
    // -x*log1pmx((M*fg-x)/x) + x*log(fg) + M - M*fg =
    //
    // Note, that fg has at most 10 bits.  If M and x are suitably
    // "nice" -- such as being integers or half-integers -- then
    // we can compute M*fg as well as x * bd0_scale[.][.] without
    // rounding errors.
    //

    add1(-x * log1pmx((m * fg - x) / x), &mut yh, &mut yl);
    if fg == 1.0 {
        return (yl, yh);
    }
    for j in 0..4 {
        add1(x * BD0_SCALE[i as usize][j] as f64, &mut yh, &mut yl);
        add1(-x * BD0_SCALE[0][j] as f64, &mut yh, &mut yl);
    }
    if !r_finite(yh) {
        yh = ML_POSINF;
        yl = 0.0;
        return (yh, yl);
    }
    add1(m, &mut yh, &mut yl);
    add1(-m * fg, &mut yh, &mut yl);
    (yh, yl)
}
