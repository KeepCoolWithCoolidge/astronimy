import angle
import coords
import interpol
import math

type
  TransitBody = enum
    StarOrPlanet
    Sun
    Moon
  TransitType = enum
    Rise
    Transit
    Set

proc m(
  transitType: TransitType,
  H0: float64,
  asc: float64,
  L: float64,
  theta0: float64
): float64 {.inline.} =
  var m = (asc + L - theta0) / angle.TWO_PI
  let p = H0 / angle.TWO_PI

  m = case transitType
      of Transit: m
      of Rise: m - p
      of Set: m + p
  
  if m < 0:
    m += 1.0
  elif m > 1.0:
    m -= 1.0
  
  m

proc time*(
  transitType: TransitType,
  transitBody: TransitBody,
  geographPoint: coords.GeographPoint,
  eqPoint1: coords.EqPoint,
  eqPoint2: coords.EqPoint,
  eqPoint3: coords.EqPoint,
  apprntGreenwhichSidr: float64,
  deltaT: float64,
  moonEqHzParallax: float64
): tuple[hours, min: int64,  sec: float64] =
  ## Computes the time of transit for a celestial body
  let h0 = case transitBody
           of StarOrPlanet: -0.5667.degToRad()
           of Sun: -0.8333.degToRad()
           of Moon: 0.7275 * moonEqHzParallax -
                    0.5667.degToRad()
  var H0 = (
    (h0.sin() - geographPoint.lat.sin() * eqPoint2.dec.sin()) /
    (geographPoint.lat.cos() * eqPoint2.dec.cos())
  ).arccos()
  H0 = angle.limitToTwoPI(H0)

  var m = m(
    transitType,
    H0,
    eqPoint2.asc,
    geographPoint.long,
    apprntGreenwhichSidr
  )
  let theta0 = apprntGreenwhichSidr + m * 360.985647.degToRad()

  let d = m + deltaT / 86400.0

  let asc = interpol.threeValues(eqPoint1.asc, eqPoint2.asc, eqPoint3.asc, d)

  let dec = case transitType
            of Transit: 0.0
            of Rise: interpol.threeValues(
              eqPoint1.dec, eqPoint2.dec,
              eqPoint3.dec, d
            )
            of Set: interpol.threeValues(
              eqPoint1.dec, eqPoint2.dec,
              eqPoint3.dec, d
            )

  var H = coords.hrAnglFrmObserverLong(
    theta0, geographPoint.long, asc
  ).radToDeg()
  H = angle.limitTo360(H)
  if H > 180.0:
    H -= 360.0
  H = H.degToRad()

  var h = case transitType
          of Transit: 0.0
          of Rise, Set: coords.altFrmEq(
            H, dec, geographPoint.lat
          )
  
  m = case transitType
      of Transit: m - H / angle.TWO_PI
      of Rise, Set: m - (h - h0) / (angle.TWO_PI * dec.cos() * geographPoint.lat.cos() * H.sin())
  
  h = 24.0 * m
  let hour = h.int64
  m = (h - (hour.float64)) * 60.0
  let minute = m.int64
  let second = (m - (minute.float64)) * 60.0

  (hour, minute, second)