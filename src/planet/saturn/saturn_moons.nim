import planet
import precess
import time
import math

type
  Moon* = enum
    Mimas
    Enceladus
    Thethys
    Dione
    Rhea
    Titan
    Hyperion
    Iapetus

type
  Info* = object
    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11: float64
    W0, W1, W2, W3, W4, W5, W6, W7, W8: float64
    s1, c1, s2, c2: float64
    e1: float64
    lambda0, beta0, delta: float64

proc createInfoStruct(JD: float64): Info =
  let angle1 = 28.0817.degToRad()
  let angle2 = 168.8112.degToRad()

  var info = Info(
    t1: JD - 2411093.0,
    t2: 0.0,
    t3: (JD - 2433282.423)/365.25 + 1950.0,
    t4: JD - 2411368.0,
    t5: 0.0,
    t6: JD - 2415020.0,
    t7: 0.0, t8: 0.0,
    t9: (JD - 2442000.5)/365.25,
    t10: JD - 2409786.0,
    t11: 0.0,
    W0: 0.0, W1: 0.0, W2: 0.0, W3: 0.0, W4: 0.0,
    W5: 0.0, W6: 0.0, W7: 0.0, W8: 0.0,
    s1: angle1.sin(),
    c1: angle1.cos(),
    s2: angle2.sin(),
    c2: angle2.cos(),
    e1: 0.0, lambda0: 0.0, beta0: 0.0, delta: 0.0
  )

  info.t2 = info.t1/365.25;
  info.t5 = info.t4/365.25;
  info.t7 = info.t6/36525.0;
  info.t8 = info.t6/365.25;
  info.t11 = info.t10/36525.0;

  info.W0 = 5.095 * (info.t3 - 1866.39).degToRad();
  info.W1 = (74.4     + 32.39*info.t2).degToRad();
  info.W2 = (134.3    + 92.62*info.t2).degToRad();
  info.W3 = (42.0     - 0.5118*info.t5).degToRad();
  info.W4 = (276.59   + 0.5118*info.t5).degToRad();
  info.W5 = (267.2635 + 1222.1136*info.t7).degToRad();
  info.W6 = (175.4762 + 1221.5515*info.t7).degToRad();
  info.W7 = (2.4891   + 0.002435*info.t7).degToRad();
  info.W8 = (113.35   - 0.2597*info.t7).degToRad();

  info.e1 = 0.05589 - 0.000346*info.t7;

  info

proc mimas(info: Info): tuple[lambda1, gamma1, omega1, r1: float64] {.inline.} =
  let L = (
      127.64 +
      381.994497*info.t1 -
      43.57*info.W0.sin() -
      0.72*(3.0*info.W0).sin() -
      0.02144*(5.0*info.W0).sin()
  ).degToRad();

  let p = (106.1 + 365.549*info.t2).degToRad();
  let M = L - p;
  let C = (
      2.18287*M.sin() +
      0.025988*(2.0*M).sin() +
      0.00043*(3.0*M).sin()
  ).degToRad();

  let lambda1 = L + C;
  let gamma1 = 1.563.degToRad();
  let omega1 = (54.5 - 365.072*info.t2).degToRad();
  let r1 = 3.06879/(1.0 + 0.01905*(M + C).cos());

  (lambda1, gamma1, omega1, r1)

proc enceladus(info: Info): tuple[lambda2, gamma2, omega2, r2: float64] {.inline.} =
  let L = (
      200.317 +
      262.7319002*info.t1 +
      0.25667*info.W1.sin() +
      0.20883*info.W2.sin()
  ).degToRad();

  let p = (309.107 + 123.44121*info.t2).degToRad();
  let M = L - p;
  let C = (
      0.55577*M.sin() +
      0.00168*(2.0*M).sin()
  ).degToRad();

  let lambda2 = L + C;
  let gamma2 = 0.0262.degToRad();
  let omega2 = (348.0 - 151.95*info.t2).degToRad();
  let r2 = 3.94118/(1.0 + 0.00485*(M + C).cos());

  (lambda2, gamma2, omega2, r2)

#proc tethys(info: Info): tuple[lambda3, gamma3, omega3, r3] =


proc apprntRectCoords*(JD: float64, moon: Moon): tuple[X, Y, Z: float64] =
  var info = createInfoStruct(JD - 0.04942)

  let (planetEclPoint, saturnEarthDist) = planet.geocentApprntEclCoords(planet.Planet.Saturn, JD)
  var (lambda0, beta0) = (planetEclPoint.long, planetEclPoint.lat)

  (lambda0, beta0) = precess.precessEclCoords(
    lambda0, beta0,
    JD,
    time.julianDay(
      time.Date(
        year: 1950,
        month: time.Month.Jan,
        decimalDay: 1.5,
        calType: time.CalType.Gregorian
      )
    )
  )

  info.lambda0 = lambda0
  info.beta0 = beta0
  info.delta = saturnEarthDist