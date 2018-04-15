import angle
import time
import coords
import math

proc nutation*(JD: float64): tuple[nutInLong, nutInOblq: float64] =
  ## Computes nutation in ecliptic longitude and obliquity
  type
    terms = object
      e1: int8
      e2: int8
      e3: int8
      e4: int8
      e5: int8
      e6: int32
      e7: int32
      e8: int32
      e9: int16
  
  let termsForNutation = [
        terms(e1:  0, e2:  0, e3:  0, e4: 0, e5: 1, e6: -171996, e7: -1742, e8: 92025,  e9: 89),
        terms(e1: -2, e2:  0, e3:  0, e4: 2, e5: 2, e6:  -13187, e7:   -16, e8:  5736,  e9: -31),
        terms(e1:  0, e2:  0, e3:  0, e4: 2, e5: 2, e6:   -2274, e7:    -2, e8:   977,  e9: -5),
        terms(e1:  0, e2:  0, e3:  0, e4: 0, e5: 2, e6:    2062, e7:     2, e8:  -895,  e9:  5),
        terms(e1:  0, e2:  1, e3:  0, e4: 0, e5: 0, e6:    1426, e7:   -34, e8:    54,  e9: -1),
        terms(e1:  0, e2:  0, e3:  1, e4: 0, e5: 0, e6:     712, e7:     1, e8:    -7,  e9:  0),
        terms(e1: -2, e2:  1, e3:  0, e4: 2, e5: 2, e6:    -517, e7:    12, e8:   224,  e9: -6),
        terms(e1:  0, e2:  0, e3:  0, e4: 2, e5: 1, e6:    -386, e7:    -4, e8:   200,  e9:  0),
        terms(e1:  0, e2:  0, e3:  1, e4: 2, e5: 2, e6:    -301, e7:     0, e8:   129,  e9: -1),
        terms(e1: -2, e2: -1, e3:  0, e4: 2, e5: 2, e6:     217, e7:    -5, e8:   -95,  e9:  3),
        terms(e1: -2, e2:  0, e3:  1, e4: 0, e5: 0, e6:    -158, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  0, e3:  0, e4: 2, e5: 1, e6:     129, e7:     1, e8:   -70,  e9:  0),
        terms(e1:  0, e2:  0, e3: -1, e4: 2, e5: 2, e6:     123, e7:     0, e8:   -53,  e9:  0),
        terms(e1:  2, e2:  0, e3:  0, e4: 0, e5: 0, e6:      63, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3:  1, e4: 0, e5: 1, e6:      63, e7:     1, e8:   -33,  e9:  0),
        terms(e1:  2, e2:  0, e3: -1, e4: 2, e5: 2, e6:     -59, e7:     0, e8:    26,  e9:  0),
        terms(e1:  0, e2:  0, e3: -1, e4: 0, e5: 1, e6:     -58, e7:    -1, e8:    32,  e9:  0),
        terms(e1:  0, e2:  0, e3:  1, e4: 2, e5: 1, e6:     -51, e7:     0, e8:    27,  e9:  0),
        terms(e1: -2, e2:  0, e3:  2, e4: 0, e5: 0, e6:      48, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3: -2, e4: 2, e5: 1, e6:      46, e7:     0, e8:   -24,  e9:  0),
        terms(e1:  2, e2:  0, e3:  0, e4: 2, e5: 2, e6:     -38, e7:     0, e8:    16,  e9:  0),
        terms(e1:  0, e2:  0, e3:  2, e4: 2, e5: 2, e6:     -31, e7:     0, e8:    13,  e9:  0),
        terms(e1:  0, e2:  0, e3:  2, e4: 0, e5: 0, e6:      29, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  0, e3:  1, e4: 2, e5: 2, e6:      29, e7:     0, e8:   -12,  e9:  0),
        terms(e1:  0, e2:  0, e3:  0, e4: 2, e5: 0, e6:      26, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  0, e3:  0, e4: 2, e5: 0, e6:     -22, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3: -1, e4: 2, e5: 1, e6:      21, e7:     0, e8:   -10,  e9:  0),
        terms(e1:  0, e2:  2, e3:  0, e4: 0, e5: 0, e6:      17, e7:    -1, e8:     0,  e9:  0),
        terms(e1:  2, e2:  0, e3: -1, e4: 0, e5: 1, e6:      16, e7:     0, e8:    -8,  e9:  0),
        terms(e1: -2, e2:  2, e3:  0, e4: 2, e5: 2, e6:     -16, e7:     1, e8:     7,  e9:  0),
        terms(e1:  0, e2:  1, e3:  0, e4: 0, e5: 1, e6:     -15, e7:     0, e8:     9,  e9:  0),
        terms(e1: -2, e2:  0, e3:  1, e4: 0, e5: 1, e6:     -13, e7:     0, e8:     7,  e9:  0),
        terms(e1:  0, e2: -1, e3:  0, e4: 0, e5: 1, e6:     -12, e7:     0, e8:     6,  e9:  0),
        terms(e1:  0, e2:  0, e3:  2, e4: -2, e5: 0, e6:      11, e7:     0, e8:     0,  e9:  0),
        terms(e1:  2, e2:  0, e3: -1, e4: 2, e5: 1, e6:     -10, e7:     0, e8:     5,  e9:  0),
        terms(e1:  2, e2:  0, e3:  1, e4: 2, e5: 2, e6:      -8, e7:     0, e8:     3,  e9:  0),
        terms(e1:  0, e2:  1, e3:  0, e4: 2, e5: 2, e6:       7, e7:     0, e8:    -3,  e9:  0),
        terms(e1: -2, e2:  1, e3:  1, e4: 0, e5: 0, e6:      -7, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2: -1, e3:  0, e4: 2, e5: 2, e6:      -7, e7:     0, e8:     3,  e9:  0),
        terms(e1:  2, e2:  0, e3:  0, e4: 2, e5: 1, e6:      -7, e7:     0, e8:     3,  e9:  0),
        terms(e1:  2, e2:  0, e3:  1, e4: 0, e5: 0, e6:       6, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  0, e3:  2, e4: 2, e5: 2, e6:       6, e7:     0, e8:    -3,  e9:  0),
        terms(e1: -2, e2:  0, e3:  1, e4: 2, e5: 1, e6:       6, e7:     0, e8:    -3,  e9:  0),
        terms(e1:  2, e2:  0, e3: -2, e4: 0, e5: 1, e6:      -6, e7:     0, e8:     3,  e9:  0),
        terms(e1:  2, e2:  0, e3:  0, e4: 0, e5: 1, e6:      -6, e7:     0, e8:     3,  e9:  0),
        terms(e1:  0, e2: -1, e3:  1, e4: 0, e5: 0, e6:       5, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2: -1, e3:  0, e4: 2, e5: 1, e6:      -5, e7:     0, e8:     3,  e9:  0),
        terms(e1: -2, e2:  0, e3:  0, e4: 0, e5: 1, e6:      -5, e7:     0, e8:     3,  e9:  0),
        terms(e1:  0, e2:  0, e3:  2, e4: 2, e5: 1, e6:      -5, e7:     0, e8:     3,  e9:  0),
        terms(e1: -2, e2:  0, e3:  2, e4: 0, e5: 1, e6:       4, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  1, e3:  0, e4: 2, e5: 1, e6:       4, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3:  1, e4: -2, e5: 0, e6:       4, e7:     0, e8:     0,  e9:  0),
        terms(e1: -1, e2:  0, e3:  1, e4: 0, e5: 0, e6:      -4, e7:     0, e8:     0,  e9:  0),
        terms(e1: -2, e2:  1, e3:  0, e4: 0, e5: 0, e6:      -4, e7:     0, e8:     0,  e9:  0),
        terms(e1:  1, e2:  0, e3:  0, e4: 0, e5: 0, e6:      -4, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3:  1, e4: 2, e5: 0, e6:       3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3: -2, e4: 2, e5: 2, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1: -1, e2: -1, e3:  1, e4: 0, e5: 0, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  1, e3:  1, e4: 0, e5: 0, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2: -1, e3:  1, e4: 2, e5: 2, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  2, e2: -1, e3: -1, e4: 2, e5: 2, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  0, e2:  0, e3:  3, e4: 2, e5: 2, e6:      -3, e7:     0, e8:     0,  e9:  0),
        terms(e1:  2, e2: -1, e3:  0, e4: 2, e5: 2, e6:      -3, e7:     0, e8:     0,  e9:  0)
  ]

  let t = time.julianCent(JD)

  let M1 = angle.limitTo360(
    134.96298 + t * (477198.867398 + t * (0.0086972 + t / 56250.0))
  ).degToRad()

  let M = angle.limitTo360(
    357.52772 + t * (35999.05034 - t * (0.0001603 + t / 300000.0))
  ).degToRad()

  let D = angle.limitTo360(
    297.85036 + t * (445267.11148 - t * (0.0019142 - t / 189474.0))
  ).degToRad()

  let F = angle.limitTo360(
    93.27191 + t * (483202.017538 - t * (0.0036825 - t / 327270.0))
  ).degToRad()

  let om = angle.limitTo360(
    125.04452 - t * (1934.136261 - t * (0.0020708 - t / 450000.0))
  ).degToRad()

  var nutInLong = 0.0
  var nutInOblq = 0.0

  var `div` = 0.0001 / 3600.0

  for x in termsForNutation:
    var arg =
      (x.e1.float64) * D +
      (x.e2.float64) * M +
      (x.e3.float64) * M1 +
      (x.e4.float64) * F +
      (x.e5.float64) * om
    
    nutInLong += ((x.e6.float64) + t * (x.e7.float64) / 10.0) * arg.sin() * `div`
    nutInOblq += ((x.e8.float64) + t * (x.e9.float64) / 10.0) * arg.cos() * `div`
  
  (nutInLong.degToRad(), nutInOblq.degToRad())

proc nutationInEqCoords*(
  eqPoint: EqPoint, 
  nutInLong: float64,
  nutInOblq: float64,
  truOblq: float64
): tuple[nutAsc, nutDec: float64] =
  ## Computes nutation in equatorial coordinates
  let (asc, dec) = (eqPoint.asc, eqPoint.dec)

  let nutAsc = nutInLong * (
    truOblq.cos() +
    truOblq.sin() * asc.sin() * dec.tan()
  ) - asc.cos() * dec.tan() * nutInOblq

  let nutDec =
    truOblq.sin() * asc.cos() * nutInLong +
    asc.sin() * nutInOblq
  
  (nutAsc, nutDec)
