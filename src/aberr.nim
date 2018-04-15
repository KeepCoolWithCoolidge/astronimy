import angle
import coords
import time
import math

proc solAberr*(R: float64): float64 {.inline.} =
  ## Computes solar aberration in ecliptic longitude
  -angle.degFrmDMS(0, 0, 20.4898).degToRad() / R

proc stellAberrInEqCoords*(stellEqPoint: EqPoint, JD: float64): tuple[abrrInAsc, abbrInDec: float64] =
  ## Computes stellar aberration in equatorial coordinates
  let t = julianCent(JD)

  let l2 = 3.1761467 + 1_012.3285546 * t
  let l3 = 1.7534703 +   628.3075849 * t
  let l4 = 6.2034809 +   334.0612431 * t
  let l5 = 0.5995465 +    52.9690965 * t
  let l6 = 0.8740168 +    21.3299095 * t
  let l7 = 5.4812939 +     7.4781599 * t
  let l8 = 5.3118863 +     3.8133036 * t
  let l1 = 3.8103444 + 8_399.6847337 * t
  let  d = 5.1984667 + 7_771.3771486 * t
  let m1 = 2.3555559 + 8_328.6914289 * t
  let  f = 1.6279052 + 8_433.4661601 * t

  var x = 0.0
  var y = 0.0
  var z = 0.0

  var A: float64
  var sinA: float64
  var cosA: float64

  A = l3
  sinA = A.sin()
  cosA = A.cos()
  x += (-1_719_914.0 - 2.0 * t) * sinA - 25.0 * cosA
  y += (25.0 - 13.0 * t) * sinA + (1_578_089.0 + 156.0 * t) * cosA
  z += (10.0 + 32.0 * t) * sinA + (684_185.0 - 358.0 * t) * cosA

  A = 2.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += (6_434.0 + 141.0 * t) * sinA + (28_007.0 - 107.0 * t) * cosA
  y += (25_697.0 - 95.0 * t) * sinA + (-5_904.0 - 130.0 * t) * cosA
  z += (11_141.0 - 48.0 * t) * sinA + (-2_559.0 -  55.0 * t) * cosA

  A = l5
  sinA = A.sin()
  cosA = A.cos()
  x += 715.0 * sinA
  y += 6.0 * sinA - 657.0 * cosA
  z += -15.0 * sinA - 282.0 * cosA

  A = l1
  sinA = A.sin()
  cosA = A.cos()
  x += 715.0 * sinA
  y += -656.0 * cosA
  z += -285.0 * cosA

  A = 3.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += (486.0 - 5.0 * t) * sinA + (-236.0 - 4.0 * t) * cosA
  y += (-216.0 - 4.0 * t) * sinA + (-446.0 + 5.0 * t) * cosA
  z += -94.0 * sinA - 193.0 * cosA

  A = f
  cosA = A.cos()
  y += 26.0 * cosA
  z += -59.0 * cosA

  A = l1 + m1
  sinA = A.sin()
  cosA = A.cos()
  x += 39.0 * sinA
  y += -36.0 * cosA
  z += -16.0 * cosA

  A = 2.0 * l5
  sinA = A.sin()
  cosA = A.cos()
  x += 33.0 * sinA - 10.0 * cosA
  y += -9.0 * sinA - 30.0 * cosA
  z += -5.0 * sinA - 13.0 * cosA

  A = 2.0 * l3 - l5
  sinA = A.sin()
  cosA = A.cos()
  x += 31.0 * sinA + cosA
  y += sinA - 28.0 * cosA
  z += -12.0 * cosA

  A = 3.0 * l3 - 8.0 * l4 + 3.0 * l5
  sinA = A.sin()
  cosA = A.cos()
  x += 8.0 * sinA - 28.0 * cosA
  y += 25.0 * sinA + 8.0 * cosA
  z += 11.0 * sinA + 3.0 * cosA

  A = 5.0 * l3 - 8.0 * l4 + 3.0 * l5
  sinA = A.sin()
  cosA = A.cos()
  x += 8.0 * sinA - 28.0 * cosA
  y += -25.0 * sinA - 8.0 * cosA
  z += -11.0 * sinA - 3.0 * cosA

  A = 2.0 * l2 - l3
  sinA = A.sin()
  cosA = A.cos()
  x += 21.0 * sinA
  y += -19.0 * cosA
  z += -8.0 * cosA

  A = l2
  sinA = A.sin()
  cosA = A.cos()
  x += -19.0 * sinA
  y += 17.0 * cosA
  z += 8.0 * cosA

  A = l7
  sinA = A.sin()
  cosA = A.cos()
  x += 17.0 * sinA
  y += -16.0 * cosA
  z += -7.0 * cosA

  A = l3 - 2.0 * l5
  sinA = A.sin()
  cosA = A.cos()
  x += 16.0 * sinA
  y += 15.0 * cosA
  z += sinA + 7.0 * cosA

  A = l8
  sinA = A.sin()
  cosA = A.cos()
  x += 16.0 * sinA
  y += sinA - 15.0 * cosA
  z += -3.0 * sinA - 6.0 * cosA

  A = l3 + l5
  sinA = A.sin()
  cosA = A.cos()
  x += 11.0 * sinA - cosA
  y += -1.0 * sinA - 10.0 * cosA
  z += -1.0 * sinA - 5.0 * cosA

  A = 2.0 * (l2 - l3)
  sinA = A.sin()
  cosA = A.cos()
  x += -11.0 * cosA
  y += -10.0 * sinA
  z += -4.0 * sinA

  A = l3 - l5
  sinA = A.sin()
  cosA = A.cos()
  x += -11.0 * sinA - 2.0 * cosA
  y += -2.0 * sinA + 9.0 * cosA
  z += -1.0 * sinA + 4.0 * cosA

  A = 4.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -7.0 * sinA - 8.0 * cosA
  y += -8.0 * sinA + 6.0 * cosA
  z += 3.0 * (cosA - sinA)

  A = 3.0 * l3 - 2.0 * l5
  sinA = A.sin()
  cosA = A.cos()
  x += -10 * sinA
  y += 9.0 * cosA
  z += 4.0 * cosA

  A = l2 - 2.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -9.0 * sinA
  y += -9.0 * cosA
  z += -4.0 * cosA

  A = 2.0 * l2 - 3.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -9.0 * sinA
  y += -8.0 * cosA
  z += -4.0 * cosA

  A = 2.0 * l6
  sinA = A.sin()
  cosA = A.cos()
  x += -9.0 * cosA
  y += -8.0 * sinA
  z += -3.0 * sinA

  A = 2.0 * l2 - 4.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -9.0 * cosA
  y += 8.0 * sinA
  z += 3.0 * sinA

  A = 3.0 * l3 - 2.0 * l4
  sinA = A.sin()
  cosA = A.cos()
  x += 8.0 * sinA
  y += -8.0 * cosA
  z += -3.0 * cosA

  A = l1 + 2.0 * d - m1
  sinA = A.sin()
  cosA = A.cos()
  x += 8.0 * sinA
  y += -7.0 * cosA
  z += -3.0 * cosA

  A = 8.0 * l2 - 12.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -4.0 * sinA - 7.0 * cosA
  y += -6.0 * sinA + 4.0 * cosA
  z += -3.0 * sinA + 2.0 * cosA

  A = 8.0 * l2 - 14.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -4.0 * sinA - 7.0 * cosA
  y += 6.0 * sinA - 4.0 * cosA
  z += 3.0 * sinA - 2.0 * cosA

  A = 2.0 * l4
  sinA = A.sin()
  cosA = A.cos()
  x += -6.0 * sinA - 5.0 * cosA
  y += -4.0 * sinA + 5.0 * cosA
  z += 2.0 * (cosA - sinA)

  A = 3.0 * l2 - 4.0 * l3
  sinA = A.sin()
  cosA = A.cos()
  x += -1.0 * (sinA + cosA)
  y += -2.0 * sinA - 7.0 * cosA
  z += sinA - 4.0 * cosA

  A = 2.0 * (l3 - l5)
  sinA = A.sin()
  cosA = A.cos()
  x += 4.0 * sinA - 6.0 * cosA
  y += -5.0 * sinA - 4.0 * cosA
  z += -2.0 * (sinA + cosA)

  A = 3.0 * (l2 - l3)
  sinA = A.sin()
  cosA = A.cos()
  x += -7.0 * cosA
  y += -6.0 * sinA
  z += -3.0 * sinA

  A = 2.0 * (l3 - l4)
  sinA = A.sin()
  cosA = A.cos()
  x += 5.0 * (sinA - cosA)
  y += -4.0 * sinA - 5.0 * cosA
  z += -2.0 * (sinA + cosA)

  A = l1 - 2.0 * d
  sinA = A.sin()
  cosA = A.cos()
  x += 5.0 * sinA
  y += -5.0 * cosA
  z += -2.0 * cosA

  let c = 17_314_463_350.0

  let (asc, dec) = (stellEqPoint.asc, stellEqPoint.dec)

  let deltaAsc = (y * asc.cos() - x * asc.sin()) / (c * dec.cos())
  let deltaDec = -(((x * asc.cos() + y * asc.sin()) * dec.sin() - z * dec.cos())) / c

  (delta_asc, delta_dec)
