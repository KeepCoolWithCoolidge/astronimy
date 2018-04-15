import angle
import coords
import time
import math
import util
import sequtils

const
  IncOfMnLunarEq*: float64 = angle.degFrmDMS(1, 32, 32.7).degToRad()
  CosIncOfMnLunarEq*: float64 = cos(IncOfMnLunarEq)
  SinIncOfMnLunarEq*: float64 = sin(IncOfMnLunarEq)

proc mnAscendNode(JC: float64): float64

proc eqHzParllx*(earthMoonDist: float64): float64 {.inline.} =
  ## Computes the equatorial horizontal parallax of the Moon
  (6378.14 / earthMoonDist).arcsin()

proc semidiameter*(earthMoonDist: float64): float64 {.inline.} =
  0.272481 * eqHzParllx(earthMoonDist).sin()

proc A(mnGeocentMoonLong : float64, 
       appGeocentMoonLat : float64, 
       longOfMnAscentNode: float64
): float64 =
  let W = mnGeocentMoonLong - longOfMnAscentNode
  (
    W.sin() * appGeocentMoonLat.cos() * CosIncOfMnLunarEq -
    appGeocentMoonLat.sin() * SinIncOfMnLunarEq
  ).arctan2(W.cos() * appGeocentMoonLat.cos())

proc F(JC: float64): float64 =
  angle.limitTo360(
    HornerEval(
      JC,
      93.272095,
      [
        483202.0175233,
            -0.0036539,
            -1.0 / 3526000.0,
             1.0 / 863310000.0
      ]
    )
  ).degToRad()

proc E(JC: float64): float64 {.inline.} =
  1.0 - JC * (0.002516 + JC * 0.0000074)

proc DMM1(JC: float64): tuple[D, M, M1: float64] {.inline.} =
  let D = angle.limitTo360(
    HornerEval(
      JC,
      297.8501921,
      [
        445267.1114034,
            -0.0018819,
             1.0 / 545868.0,
            -1.0 / 113065000.0
      ]
    )
  ).degToRad()

  let M = angle.limitTo360(
    HornerEval(
      JC,
      357.5291092,
      [
        35999.0502909,
           -0.0001536,
            1.0 / 24490000.0
      ]
    )
  ).degToRad()

  let M1 = angle.limitTo360(
    HornerEval(
      JC,
      134.9633964,
      [
        447198.8675055,
             0.0087414,
             1.0 / 69699.0,
            -1.0 / 14712000.0
      ]
    )
  ).degToRad()

  (D, M, M1)

proc rho_sig(D, M1, F: float64): tuple[rho, sig: float64] {.inline.} =
  let x2F = 2.0 * F
  let x2F2D = x2F - 2.0 * D

  (
    (
      - 0.02752 * M1.cos() -
        0.02245 * F.sin() +
        0.00684 * (M1 - x2F).cos() -
        0.00293 * (x2F).cos() -
        0.00085 * (x2F2D).cos() -
        0.00054 * (M1 - 2.0 * D).cos() -
        0.00020 * (
          (M1 + F).sin() +
          (M1 + x2F).cos() +
          (M1 - F).cos()
        ) +
      0.00014 * (M1 + x2F2D).cos()
    ).degToRad(),

    (
      - 0.02816 * M1.sin() +
        0.02244 * F.cos() -
        0.00682 * (M1 - x2F).sin() -
        0.00279 * (x2F).sin() -
        0.00083 * (x2F2D).sin() +
        0.00069 * (M1 - 2.0 * D).sin() +
        0.0004  * (M1 + F).cos() -
        0.00025 * (2.0 * M1).sin() -
        0.00023 * (M1 + x2F).sin() +
        0.0002  * (M1 - F).cos() +
        0.00019 * (M1 - F).sin() +
        0.00013 * (M1 + x2F2D).sin() -
        0.0001  * (M1 - 3.0 * F).cos()
    ).degToRad()
  )

proc opticalLibr*(
  JD              : float64,
  mnEclLongMoon   : float64,
  apprntEclLatMoon: float64
): tuple[opticalLibrInLong, opticalLibrInLat: float64] =
  ## Computes the optical librations of the Moon in longitude and
  ## latitude

  let JC = time.julianCent(JD)
  let F = F(JC)

  let longOfMnAscendNode = mnAscendNode(JC)
  let W = mnEclLongMoon - longOfMnAscendNode

  let A = A(
    mnEclLongMoon, apprntEclLatMoon, longOfMnAscendNode
  )

  let b1 = (
    - W.sin() * apprntEclLatMoon.cos() * SinIncOfMnLunarEq - 
    apprntEclLatMoon.sin() * CosIncOfMnLunarEq
  ).arcsin()

  (angle.limitToTwoPI(A - F), b1)

proc physicalLibr*(
  JD              : float64,
  mnEclLongMoon   : float64,
  apprntEclLatMoon: float64,
  opticalLibLat   : float64
): tuple[physicalLibrInLong, physicalLibrInLat: float64] =
  ## Computes the physical librations of the Moon in longitude and
  ## latitude
  let JC = time.julianCent(JD)
  let K1 = (119.75 + 131.849 * JC).degToRad()
  let K2 = (72.560 +  20.186 * JC).degToRad()

  let longOfMnAscendNode = mnAscendNode(JC)
  let (D, M, M1) = DMM1(JC)
  let F = F(JC)
  let E = E(JC)

  let x2F = 2.0 * F
  let x2D = 2.0 * D
  let x2F2D = x2F - x2D

  let (rho, sig) = rho_sig(D, M1, F)
  let tau = (
    0.02520 * E * M.sin() +
    0.00473 * (2.0 * M1 - x2F).sin() -
    0.00467 * M1.sin() +
    0.00396 * K1.sin() +
    0.00276 * (2.0 * M1 - x2D).sin() +
    0.00196 * longOfMnAscendNode.sin() -
    0.00183 * (M1 - F).cos() +
    0.00115 * (M1 - x2D).sin() -
    0.00096 * (M1 - D).sin() +
    0.00046 * (x2F2D).sin() -
    0.00039 * (M1 - F).sin() -
    0.00032 * (M1 - M - D).sin() +
    0.00027 * (2.0 * M1 - M - x2D).sin() +
    0.00023 * K2.sin() -
    0.00014 * (
      x2D.sin() -
      (2.0 * M1 - x2F).cos()
    ) -
    0.00012 * (
      (M1 - x2F).sin() +
      (2.0 * M1).sin()
    ) +
    0.0001 * (
      2.0 * (M1 - M - D)).sin()
  ).degToRad()

  let A = A(
    mnEclLongMoon, apprntEclLatMoon, longOfMnAscendNode
  )

  (
    opticalLibLat.tan() * (rho * A.cos() + sig*A.sin()) - tau,
    sig * A.cos() - rho * A.sin()
  )

proc totalLibr*(
  JD              : float64,
  mnEclLongMoon   : float64,
  apprntEclLatMoon: float64
): tuple[totalLibrInLong, totalLibrInLat: float64] =
  ## Computes the total librations of the Moon in longitude and
  ## latitude
  let (optLong, optLat) = opticalLibr(
    mnEclLongMoon, apprntEclLatMoon, JD
  )
  let (physLong, physLat) = physicalLibr(
    mnEclLongMoon, apprntEclLatMoon, JD, optLat
  )

  (optLong + physLong, optLat + physLat)

proc posAngleOfAxisOfRot*(
  JD              : float64,
  mnAscentNodeLong: float64,
  totalLibLat     : float64,
  nutInLong       : float64,
  trueOblq        : float64,
  apprntMoonAsc   : float64
): float64 =
  ## Computes the position angle of the axis of rotation of the Moon
  let JC = time.julianCent(JD)
  let (D, _, M1) = DMM1(JC)
  let F = F(JC)
  let (rho, sig) = rho_sig(D, M1, F)

  let V = mnAscentNodeLong + nutInLong + sig / SinIncOfMnLunarEq
  let X = (IncOFMnLunarEq + rho).sin() * V.sin()
  let Y =
    (IncOfMnLunarEq + rho).sin() * V.cos() * trueOblq.cos() -
    (IncOfMnLunarEq + rho).cos() * trueOblq.sin()
  let w = X.arctan2(Y)

  (
    (X * X + Y * Y).sqrt() * (apprntMoonAsc - w).cos() / totalLibLat.cos()
  ).arcsin()

proc topocentLibrByDiffCorrections*(
  observerLat        : float64,
  geocentDecMoon     : float64,
  localHourAngl      : float64,
  geocentHzParllxMoon: float64,
  posAnglAxisOfRot   : float64,
  totalLibLat        : float64
): tuple[topocentLibrInLong, topocentLibrInLat, topcentLibrInP: float64] =
  ## Computes the topocentric librations of the Moon
  let Q = (observerLat.cos() * localHourAngl.sin()).arctan2(
    geocentDecMoon.cos() * observerLat.sin() -
    geocentDecMoon.sin() * observerLat.cos() * localHourAngl.cos()
  )
  let z = (
    geocentDecMoon.sin() * observerLat.sin() +
    geocentDecMoon.cos() * observerLat.cos() * localHourAngl.cos()
  ).arccos()
  let pi1 = geocentHzParllxMoon * (z.sin() + 0.0084 * (2.0 * z).sin())

  let deltaL = -pi1 * (Q - posAnglAxisOfRot).sin() / totalLibLat.cos()
  let deltaB = -pi1 * (Q - posAnglAxisOfRot).cos()
  let deltaP = 
    deltaL * (totalLibLat + deltaB).sin() -
    pi1 * Q.sin() * geocentDecMoon.tan()
  
  (deltaL, deltaB, deltaP)

proc geocentEclPos*(JD: float64): tuple[moonEclPoint: EclPoint, radVec: float64] =
  ## Computes the geocentric ecliptic position of the Moon,
  ## referred to the mean equinox of the date
  let JC = time.julianCent(JD)
  let (D, M, M1) = DMM1(JC)
  let F = F(JC)
  let E = E(JC)
  let L1 = angle.limitTo360(
    HornerEval(
      JC,
      218.3164477,
      [
        481267.88123421,
            -0.0015786,
             1.0 / 538841.0,
            -1.0 / 65194000.0
      ]
    )
  ).degToRad()

  let A1 = angle.limitTo360(119.75 + 131.849 * JC).degToRad()
  let A2 = angle.limitTo360(53.090 + 479264.29 * JC).degToRad()
  let A3 = angle.limitTo360(313.45 + 481266 * JC).degToRad()

  type
    lrterms = object
      e1: int8
      e2: int8
      e3: int8
      e4: int8
      e5: int32
      e6: int32
    bterms = object
      e1: int8
      e2: int8
      e3: int8
      e4: int8
      e5: int32
  
  let termsForLr = [
    lrterms(e1: 0, e2:  0, e3:  1, e4:  0, e5:  6288774, e6: -20905355),
    lrterms(e1: 2, e2:  0, e3: -1, e4:  0, e5:  1274027, e6: -3699111),
    lrterms(e1: 2, e2:  0, e3:  0, e4:  0, e5:  658314,  e6: -2955968),
    lrterms(e1: 0, e2:  0, e3:  2, e4:  0, e5:  213618,  e6: -569925),
    lrterms(e1: 0, e2:  1, e3:  0, e4:  0, e5: -185116,  e6:  48888),
    lrterms(e1: 0, e2:  0, e3:  0, e4:  2, e5: -114332,  e6: -3149),
    lrterms(e1: 2, e2:  0, e3: -2, e4:  0, e5:  58793,   e6:  246158),
    lrterms(e1: 2, e2: -1, e3: -1, e4:  0, e5:  57066,   e6: -152138),
    lrterms(e1: 2, e2:  0, e3:  1, e4:  0, e5:  53322,   e6: -170733),
    lrterms(e1: 2, e2: -1, e3:  0, e4:  0, e5:  45758,   e6: -204586),
    lrterms(e1: 0, e2:  1, e3: -1, e4:  0, e5: -40923,   e6: -129620),
    lrterms(e1: 1, e2:  0, e3:  0, e4:  0, e5: -34720,   e6:  108743),
    lrterms(e1: 0, e2:  1, e3:  1, e4:  0, e5: -30383,   e6:  104755),
    lrterms(e1: 2, e2:  0, e3:  0, e4: -2, e5:  15327,   e6:  10321),
    lrterms(e1: 0, e2:  0, e3:  1, e4:  2, e5: -12528,   e6:  0),
    lrterms(e1: 0, e2:  0, e3:  1, e4: -2, e5:  10980,   e6:  79661),
    lrterms(e1: 4, e2:  0, e3: -1, e4:  0, e5:  10675,   e6: -34782),
    lrterms(e1: 0, e2:  0, e3:  3, e4:  0, e5:  10034,   e6: -23210),
    lrterms(e1: 4, e2:  0, e3: -2, e4:  0, e5:  8548,    e6: -21636),
    lrterms(e1: 2, e2:  1, e3: -1, e4:  0, e5: -7888,    e6:  24208),
    lrterms(e1: 2, e2:  1, e3:  0, e4:  0, e5: -6766,    e6:  30824),
    lrterms(e1: 1, e2:  0, e3: -1, e4:  0, e5: -5163,    e6: -8379),
    lrterms(e1: 1, e2:  1, e3:  0, e4:  0, e5:  4987,    e6: -16675),
    lrterms(e1: 2, e2: -1, e3:  1, e4:  0, e5:  4036,    e6: -12831),
    lrterms(e1: 2, e2:  0, e3:  2, e4:  0, e5:  3994,    e6: -10445),
    lrterms(e1: 4, e2:  0, e3:  0, e4:  0, e5:  3861,    e6: -11650),
    lrterms(e1: 2, e2:  0, e3: -3, e4:  0, e5:  3665,    e6:  14403),
    lrterms(e1: 0, e2:  1, e3: -2, e4:  0, e5: -2689,    e6: -7003),
    lrterms(e1: 2, e2:  0, e3: -1, e4:  2, e5: -2602,    e6:  0),
    lrterms(e1: 2, e2: -1, e3: -2, e4:  0, e5:  2390,    e6:  10056),
    lrterms(e1: 1, e2:  0, e3:  1, e4:  0, e5: -2348,    e6:  6322),
    lrterms(e1: 2, e2: -2, e3:  0, e4:  0, e5:  2236,    e6: -9884),
    lrterms(e1: 0, e2:  1, e3:  2, e4:  0, e5: -2120,   e6:  5751),
    lrterms(e1: 0, e2:  2, e3:  0, e4:  0, e5: -2069,   e6:  0),
    lrterms(e1: 2, e2: -2, e3: -1, e4:  0, e5:  2048,   e6: -4950),
    lrterms(e1: 2, e2:  0, e3:  1, e4: -2, e5: -1773,   e6:  4130),
    lrterms(e1: 2, e2:  0, e3:  0, e4:  2, e5: -1595,   e6:  0),
    lrterms(e1: 4, e2: -1, e3: -1, e4:  0, e5:  1215,   e6: -3958),
    lrterms(e1: 0, e2:  0, e3:  2, e4:  2, e5: -1110,   e6:  0),
    lrterms(e1: 3, e2:  0, e3: -1, e4:  0, e5: -892,     e6:  3258),
    lrterms(e1: 2, e2:  1, e3:  1, e4:  0, e5: -810,     e6:  2616),
    lrterms(e1: 4, e2: -1, e3: -2, e4:  0, e5:  759,     e6: -1897),
    lrterms(e1: 0, e2:  2, e3: -1, e4:  0, e5: -713,     e6: -2117),
    lrterms(e1: 2, e2:  2, e3: -1, e4:  0, e5: -700,     e6:  2354),
    lrterms(e1: 2, e2:  1, e3: -2, e4:  0, e5:  691,     e6:  0),
    lrterms(e1: 2, e2: -1, e3:  0, e4: -2, e5:  596,     e6:  0),
    lrterms(e1: 4, e2:  0, e3:  1, e4:  0, e5:  549,     e6: -1423),
    lrterms(e1: 0, e2:  0, e3:  4, e4:  0, e5:  537,     e6: -1117),
    lrterms(e1: 4, e2: -1, e3:  0, e4:  0, e5:  520,     e6: -1571),
    lrterms(e1: 1, e2:  0, e3: -2, e4:  0, e5: -487,     e6: -1739),
    lrterms(e1: 2, e2:  1, e3:  0, e4: -2, e5: -399,     e6:  0),
    lrterms(e1: 0, e2:  0, e3:  2, e4: -2, e5: -381,     e6: -4421),
    lrterms(e1: 1, e2:  1, e3:  1, e4:  0, e5:  351,     e6:  0),
    lrterms(e1: 3, e2:  0, e3: -2, e4:  0, e5: -340,     e6:  0),
    lrterms(e1: 4, e2:  0, e3: -3, e4:  0, e5:  330,     e6:  0),
    lrterms(e1: 2, e2: -1, e3:  2, e4:  0, e5:  327,     e6:  0),
    lrterms(e1: 0, e2:  2, e3:  1, e4:  0, e5: -323,     e6:  1165),
    lrterms(e1: 1, e2:  1, e3: -1, e4:  0, e5:  299,     e6:  0),
    lrterms(e1: 2, e2:  0, e3:  3, e4:  0, e5:  294,     e6:  0),
    lrterms(e1: 2, e2:  0, e3: -1, e4: -2, e5:  0,       e6:  8752)
  ]

  let termsForB = [
    bterms(e1: 0, e2:  0, e3:  0, e4:  1, e5: 5128122),
    bterms(e1: 0, e2:  0, e3:  1, e4:  1, e5:  280602),
    bterms(e1: 0, e2:  0, e3:  1, e4: -1, e5:  277693),
    bterms(e1: 2, e2:  0, e3:  0, e4: -1, e5:  173237),
    bterms(e1: 2, e2:  0, e3: -1, e4:  1, e5:  55413),
    bterms(e1: 2, e2:  0, e3: -1, e4: -1, e5:  46271),
    bterms(e1: 2, e2:  0, e3:  0, e4:  1, e5:  32573),
    bterms(e1: 0, e2:  0, e3:  2, e4:  1, e5:  17198),
    bterms(e1: 2, e2:  0, e3:  1, e4: -1, e5:  9266),
    bterms(e1: 0, e2:  0, e3:  2, e4: -1, e5:  8822),
    bterms(e1: 2, e2: -1, e3:  0, e4: -1, e5:  8216),
    bterms(e1: 2, e2:  0, e3: -2, e4: -1, e5:  4324),
    bterms(e1: 2, e2:  0, e3:  1, e4:  1, e5:  4200),
    bterms(e1: 2, e2:  1, e3:  0, e4: -1, e5: -3359),
    bterms(e1: 2, e2: -1, e3: -1, e4:  1, e5:  2463),
    bterms(e1: 2, e2: -1, e3:  0, e4:  1, e5:  2211),
    bterms(e1: 2, e2: -1, e3: -1, e4: -1, e5:  2065),
    bterms(e1: 0, e2:  1, e3: -1, e4: -1, e5: -1870),
    bterms(e1: 4, e2:  0, e3: -1, e4: -1, e5:  1828),
    bterms(e1: 0, e2:  1, e3:  0, e4:  1, e5: -1794),
    bterms(e1: 0, e2:  0, e3:  0, e4:  3, e5: -1749),
    bterms(e1: 0, e2:  1, e3: -1, e4:  1, e5: -1565),
    bterms(e1: 1, e2:  0, e3:  0, e4:  1, e5: -1491),
    bterms(e1: 0, e2:  1, e3:  1, e4:  1, e5: -1475),
    bterms(e1: 0, e2:  1, e3:  1, e4: -1, e5: -1410),
    bterms(e1: 0, e2:  1, e3:  0, e4: -1, e5: -1344),
    bterms(e1: 1, e2:  0, e3:  0, e4: -1, e5: -1335),
    bterms(e1: 0, e2:  0, e3:  3, e4:  1, e5:  1107),
    bterms(e1: 4, e2:  0, e3:  0, e4: -1, e5:  1021),
    bterms(e1: 4, e2:  0, e3: -1, e4:  1, e5:  833),
    bterms(e1: 0, e2:  0, e3:  1, e4: -3, e5:  777),
    bterms(e1: 4, e2:  0, e3: -2, e4:  1, e5:  671),
    bterms(e1: 2, e2:  0, e3:  0, e4: -3, e5:  607),
    bterms(e1: 2, e2:  0, e3:  2, e4: -1, e5:  596),
    bterms(e1: 2, e2: -1, e3:  1, e4: -1, e5:  491),
    bterms(e1: 2, e2:  0, e3: -2, e4:  1, e5: -451),
    bterms(e1: 0, e2:  0, e3:  3, e4: -1, e5:  439),
    bterms(e1: 2, e2:  0, e3:  2, e4:  1, e5:  422),
    bterms(e1: 2, e2:  0, e3: -3, e4: -1, e5:  421),
    bterms(e1: 2, e2:  1, e3: -1, e4:  1, e5: -366),
    bterms(e1: 2, e2:  1, e3:  0, e4:  1, e5: -351),
    bterms(e1: 4, e2:  0, e3:  0, e4:  1, e5:  331),
    bterms(e1: 2, e2: -1, e3:  1, e4:  1, e5:  315),
    bterms(e1: 2, e2: -2, e3:  0, e4: -1, e5:  302),
    bterms(e1: 0, e2:  0, e3:  1, e4:  3, e5: -283),
    bterms(e1: 2, e2:  1, e3:  1, e4: -1, e5: -229),
    bterms(e1: 1, e2:  1, e3:  0, e4: -1, e5:  223),
    bterms(e1: 1, e2:  1, e3:  0, e4:  1, e5:  223),
    bterms(e1: 0, e2:  1, e3: -2, e4: -1, e5: -220),
    bterms(e1: 2, e2:  1, e3: -1, e4: -1, e5: -220),
    bterms(e1: 1, e2:  0, e3:  1, e4:  1, e5: -185),
    bterms(e1: 2, e2: -1, e3: -2, e4: -1, e5:  181),
    bterms(e1: 0, e2:  1, e3:  2, e4:  1, e5: -177),
    bterms(e1: 4, e2:  0, e3: -2, e4: -1, e5:  176),
    bterms(e1: 4, e2: -1, e3: -1, e4: -1, e5:  166),
    bterms(e1: 1, e2:  0, e3:  1, e4: -1, e5: -164),
    bterms(e1: 4, e2:  0, e3:  1, e4: -1, e5:  132),
    bterms(e1: 1, e2:  0, e3: -1, e4: -1, e5: -119),
    bterms(e1: 4, e2: -1, e3:  0, e4: -1, e5:  115),
    bterms(e1: 2, e2: -2, e3:  0, e4:  1, e5:  107)
  ]

  var l = 0.0
  var r = 0.0
  var b = 0.0

  for x in termsForLr:
    var arg =
      x.e1.float64 * D +
      x.e2.float64 * M +
      x.e3.float64 * M1 +
      x.e4.float64 * F
    
    var t = if (x.e1).abs() == 1: E
            elif (x.e1).abs() == 2: E * E
            else: 1.0
    l += x.e5.float64 * t * arg.sin()
    r += x.e6.float64 * t * arg.cos()

  for x in termsForB:
    var arg =
      x.e1.float64 * D +
      x.e2.float64 * M +
      x.e3.float64 * M1 +
      x.e4.float64 * F

    var t = x.e5.float64 * arg.sin()

    t  = if (x.e2).abs() == 1: t * E
         elif (x.e2).abs() == 2: t * E * E
         else: t
    
    b += t
  
  l += 
      3958.0 * A1.sin() +
      1962.0 * (L1 - F).sin() +
        318.0 * A2.sin()
  
  b +=
     -2235.0 * L1.sin() +
        382.0 * A3.sin() +
        175.0 * (
          (A1 - F).sin() +
          (A1 + F).sin()
        ) +
      127.0 * (L1 - M1).sin() -
      115.0 * (L1 + M1).sin()

  l = l.degToRad()
  b = b.degToRad()

  let eclPoint = coords.EclPoint(
    long: L1 + l/1000000.0,
    lat: b / 1000000.0
  )

  let radVec = 385000.56 + r / 1000.0

  (eclPoint, radVec)

proc mnAscendNode(JC: float64): float64 =
  ## Computes the longitude of the mean ascending node of the Moon
  angle.limitTo360(
    HornerEval(
      JC,
      125.0445479,
      [
        -1934.1362891,
            0.020754,
            1.0 / 467441.0,
           -1.0 / 60616000.0
      ]
    )
  ).degToRad()

proc trueAscendNode*(JC: float64): float64 =
  ## Computes the longitude of the true ascending node of the Moon
  let (D, M, M1) = DMM1(JC)
  let F = F(JC)

  mnAscendNode(JC) + (
    -1.4979 * (2.0 * (D - F)).sin() -
     0.1500 * M.sin() -
     0.1226 * (2.0 * D).sin() +
     0.1176 * (2.0 * F).sin() -
     0.0801 * (2.0 * (M1 - F)).sin()
  ).degToRad()

proc mnPerigee*(JC: float64): float64 =
  ## Computes the longitude of the mean perigee of the Moon
  angle.limitTo360(
    HornerEval(
      JC,
      83.3532465,
      [
        4069.0137287,
          -0.01032,
          -1.0 / 80053.0,
           1.0 / 18999000.0
      ]
    )
  ).degToRad()

proc brightLimb*(
  sunEqPoint: EqPoint,
  moonEqPoint: EqPoint
): float64 =
  ## Computes the position angle of the Moon's bright limb
  let a = sunEqPoint.dec.cos()
  let n = a * (sunEqPoint.asc - moonEqPoint.asc).sin()
  let d =
    sunEqPoint.dec.sin() * moonEqPoint.dec.cos() -
    moonEqPoint.dec.sin() * (
      sunEqPoint.asc - moonEqPoint.asc
    ).cos() * a
  
  n.arctan2(d)

proc illuminatedFrac(
  moonGeocentElong: float64,
  earthMoonDist: float64,
  earthSunDist: float64
): float64 {.inline.} =
 let i = (earthSunDist * moonGeocentElong.sin()).arctan2(
   earthMoonDist - earthMoonDist * moonGeocentElong.cos()
 )

 (1.0 + i.cos()) / 2.0

proc illumFracFrmEqCoords*(
  sunEqPoint: EqPoint,
  moonEqPoint: EqPoint,
  earthMoonDist: float64,
  earthSunDist: float64
): float64 =
  ## Computes the illuminated fraction of the lunar disk, using equatorial
  ## coordinates
  illuminatedFrac(
    sunEqPoint.anglrSepr(moonEqPoint).arccos(),
    earthMoonDist,
    earthSunDist
  )

proc illumFracFrmEclCoords*(
  moonLong: float64,
  moonLat: float64,
  sunLong: float64,
  earthMoonDist: float64,
  earthSunDist: float64
): float64 =
  ## Computes the illuminated fraction of the lunar disk, using ecliptic
  ## coordinates
  illuminatedFrac(
    (moonLat.cos() * (moonLong - sunLong).cos()).arccos(),
    earthMoonDist,
    earthSunDist
  )

proc timeOfPassageThroughNode(k, T: float64): float64 =
  let D = (
    183.638 +
    331.73735682 * k +
    T * T * (
      0.0014852 +
      T * (0.00000209 - T * 0.00000001)
    )
  ).degToRad()

  let DTimes2 = 2.0 * D
  let M = (17.4006 + 26.8203725 * k + T * T * (0.0001186 + T * 0.00000006)).degToRad()
  let M1 = (38.3776 + 355.52747313 * k + T * T * (0.0123499 + T * (0.000014627 - T * 0.000000069))).degToRad()
  let sigma = (123.9767 - 1.44098956 * k + T * T * (0.0020608 + T * (0.0000214 - T * 0.000000016))).degToRad()
  let P = sigma + (272.75 - T * 2.3).degToRad()
  let V = (299.75 + T * (132.85 - T * 0.009173)).degToRad()

  (
    2451565.16 +
    27.212220817 * k + 
    T * T * (
      0.0002762 +
      T * (0.000000021 - T * 0.000000000088)
    ) -
    0.4721 * M1.sin() -
    0.1649 * DTimes2.sin() -
    0.0868 * (DTimes2 - M1).sin() +
    0.0084 * (DTimes2 + M1).sin() -
    0.0083 * (DTimes2 - M).sin() -
    0.0039 * (DTimes2 - M - M1).sin() +
    0.0034 * (2.0 * M1).sin() -
    0.0031 * (2.0 * (D - M1)).sin() +
    0.0030 * (DTimes2 + M).sin() +
    0.0028 * (M - M1).sin() +
    0.0026 * M.sin() +
    0.0025 * (2.0 * DTimes2).sin() +
    0.0024 * D.sin() +
    0.0022 * (M + M1).sin() +
    0.0017 * sigma.sin() +
    0.0014 * (2.0 * DTimes2 - M1).sin() +
    0.0005 * (DTimes2 + M - M1).sin() +
    0.0004 * (DTimes2 - M + M1).sin() +
    0.0003 * (
      (2.0 * (M - D)).sin() +
      (2.0 * DTimes2 - M).sin() +
      V.sin() +
      P.sin()
    )
  )

proc timeOfPassageThroughNodes*(date: Date): 
  tuple[timeOfAscendNode, timeOfDescendNode: float64] {.inline.} =
  ## Computes the times of passage of the Moon through the ascending and
  ## descending nodes, close to a given date
  let k = 13.4223 *  (decimalYear(date) - 2000.05)
  let T = k / 1342.23
  let k1 = (k.int32).float64
  let k2 = (k1.float64) + 0.5
  
  (timeOfPassageThroughNode(k1, T), timeOfPassageThroughNode(k2, T))

type
  Phase = enum
    New
    First
    Full
    Last

proc timeToPhase*(date: Date, phase: Phase): float64 =
  ## Computes the Julian day corresponding to one of the four phases
  ## of the Moon
  var K = 12.3685 * (time.decimalYear(date) - 2000.0)
  K = (K.int64).float64

  let k = case phase
          of New: 
            K
          of First:
            K + 0.25
          of Full:
            K + 0.5
          of Last:
            K + 0.75

  let T = k / 1236.85
  
  var JD = 
    2451550.09766 +
    k * 29.530588861 +
    T * HornerEval(
      T,
      0.0,
      [
        0.00015437,
       -0.00000015,
        0.00000000073
      ]
    )
  
  let E = E(T)

  let M = (
    2.5534 +
    k * 29.1053567 +
    T * HornerEval(
      T,
      0.0,
      [
        -0.0000014,
        -0.00000011
      ]
    )
  ).degToRad()

  let M1 = (
    201.5643 +
    k * 385.81693528 +
    T * HornerEval(
      T,
      0.0,
      [
        0.0107582,
        0.00001238,
       -0.000000058
      ]
    )
  ).degToRad()

  let F = (
    160.7108 +
    k * 390.67050284 +
    T * HornerEval(
      T,
      0.0,
      [
        -0.0016118,
        -0.00000227,
        0.000000011
      ]
    )
  ).degToRad()

  let omega = (
    124.7746 -
    k * 1.56375588 +
    T * HornerEval(
      T,
      0.0,
      [
        0.0020672,
        0.00000215
      ]
    )
  ).degToRad()

  let A1 = (299.77 + 0.107408 * k - 0.009173 * T * T).degToRad()
  let A2 = (251.88 + 0.016321 * k).degToRad()
  let A3 = (251.83 + 26.651886 * k).degToRad()
  let A4 = (349.42 + 36.412478 * k).degToRad()
  let A5 = (84.66 + 18.206239 * k).degToRad()
  let A6 = (141.74 + 53.303771 * k).degToRad()
  let A7 = (207.14 + 2.453732 * k).degToRad()
  let A8 = (154.84 + 7.306860 * k).degToRad()
  let A9 = (34.52 + 27.261239 * k).degToRad()
  let A10 = (207.19 + 0.121824 * k).degToRad()
  let A11 = (291.34 + 1.844379 * k).degToRad()
  let A12 = (161.72 + 24.198154 * k).degToRad()
  let A13 = (239.56 + 25.513099 * k).degToRad()
  let A14 = (331.55 + 3.592518 * k).degToRad()

  let isQuarter = case phase
                  of Last, First:
                    true
                  else:
                    false
  
  if isQuarter:
    let W =
      0.00306 -
      0.00038 * E * M.cos() +
      0.00026 * M1.cos() -
      0.00002 * ((M1 - M).cos() - (M1 + M).cos() - (2.0 * F).cos())
    
    JD  = case phase
          of Last: JD - W
          of First: JD + W
          else: JD
    
    let corrections = [
      [-0.62801,         M1],
      [ 0.17172 * E,     M],
      [-0.01183 * E,     M1 + M],
      [ 0.00862,         2.0 * M1],
      [ 0.00804,         2.0 * F],
      [ 0.00454 * E,     M1 - M],
      [ 0.00204 * E * E, 2.0 * M],
      [-0.0018,          M1 - 2.0 * F],
      [-0.0007,          M1 + 2.0 * F],
      [-0.0004,          3.0 * M1],
      [-0.00034,         2.0 * M1 - M],
      [ 0.00032 * E,     M + 2.0 * F],
      [ 0.00032 * E,     M - 2.0 * F],
      [-0.00028 * E * E, M1 + 2.0 * M],
      [ 0.00027 * E,     2.0 * M1 + M],
      [-0.00017,         omega],
      [-0.00005,         M1 - M - 2.0 * F],
      [ 0.00004,         2.0 * (M1 + F)],
      [-0.00004,         M1 + M + 2.0 * F],
      [ 0.00004,         M1 - 2.0 * M],
      [ 0.00003,         M1 + M - 2.0 * F],
      [ 0.00003,         3.0 * M],
      [ 0.00002,         2.0 * (M1 - F)],
      [ 0.00002,         M1 - M + 2.0 * F],
      [-0.00002,         3.0 * M1 + M]
    ]

    for x in corrections:
      JD += x[0] * x[1].sin()
  else:
    let isNew = case phase:
                of New: true
                else: false
    
    let sineArguments = [
      M1,
      M,
      2.0 * M1,
      2.0 * F,
      M1 - M,
      M1 + M,
      2.0 * M,
      M1 - 2.0 * F,
      M1 + 2.0 * F,
      2.0 * M1 + M,
      3.0 * M1,
      M + 2.0 * F,
      M - 2.0 * F,
      2.0 * M1 - M,
      omega,
      M1 + 2.0 * M,
      2.0 * (M1 - F),
      3.0 * M,
      M1 + M - 2.0 * F,
      2.0 * (M1 + F),
      M1 + M + 2.0 * F,
      M1 - M + 2.0 * F,
      M1 - M - 2.0 * F,
      3.0 * M1 + M,
      4.0 * M1
    ]

    let multipliers = 
      if isNew: [
        -0.4072,
         0.17241 * E,
         0.01608,
         0.01039,
         0.00739 * E,
        -0.00514 * E,
         0.00208 * E * E,
        -0.00111,
        -0.00057,
         0.00056 * E,
        -0.00042,
         0.00042 * E,
         0.00038 * E,
        -0.00024 * E,
        -0.00017,
        -0.00007,
         0.00004,
         0.00004,
         0.00003,
         0.00003,
        -0.00003,
         0.00003,
        -0.00002,
        -0.00002,
         0.00002
      ] else: [
        -0.40614,
         0.17302 * E,
         0.01614,
         0.01043,
         0.00734 * E,
        -0.00515 * E,
         0.00209 * E * E,
        -0.00111,
        -0.00057,
         0.00056 * E,
        -0.00042,
         0.00042 * E,
         0.00038 * E,
        -0.00024 * E,
        -0.00017,
        -0.00007,
         0.00004,
         0.00004,
         0.00003,
         0.00003,
        -0.00003,
         0.00003,
        -0.00002,
        -0.00002,
         0.00002
      ]

    for mx in multipliers.zip(sineArguments):
      JD += mx[0] * (mx[1]).sin()

    let additionalCorrections = [
      [0.000325, A1],
      [0.000165, A2],
      [0.000164, A3],
      [0.000126, A4],
      [0.000110, A5],
      [0.000062, A6],
      [0.00006,  A7],
      [0.000056, A8],
      [0.000047, A9],
      [0.000042, A10],
      [0.00004,  A11],
      [0.000037, A12],
      [0.000035, A13],
      [0.000023, A14]
    ]

    for x in additionalCorrections:
      JD += x[0] * x[1].sin()
    
  JD