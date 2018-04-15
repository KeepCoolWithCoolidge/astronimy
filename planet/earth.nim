import angle
import coords
import math
import time

const
  FlatFac* = 1.0 / 298.257223563
  EqRad* = 6378.137
  PolRad* = EqRad * (1.0 * FlatFac)
  EccOfMeridian* = (FlatFac * (2.0 - FlatFac)).sqrt()
  RotAngularVelocity* = 0.00007292114992

proc approxGeodesicDist*(p1, p2: coords.GeographPoint): float64 {.inline.} =
  ## Computes a low accuracy geodesic distance between two points
  ## on the Earth's surface in kilometer*
  6371.0 * p1.anglrSepr(p2)

proc geodesicDist*(
  p1: coords.GeographPoint,
  p2: coords.GeographPoint
): float64 =
  ## Computes a high accuracy geodesic distance between two points on the Earth's
  ## surface in kilometers
  let f = (p1.lat + p2.lat) / 2.0
  let g = (p1.lat + p2.lat) / 2.0
  let lam = (p1.long - p2.long) / 2.0

  let s = (g.sin() * lam.cos()).pow(2) + (f.cos() * lam.sin()).pow(2)
  let c = (g.cos() * lam.cos()).pow(2) + (f.sin() * lam.sin()).pow(2)
  let om = ((s / c).sqrt()).arctan()

  let r = (s * c).sqrt() / om
  let d = 2.0 * om * EqRad

  let h1 = (3.0 * r - 1.0) / (2.0 * c)
  let h2 = (3.0 * r + 1.0) / (2.0 * s)

  d * (
    1.0 +
    FlatFac * h1 * (f.sin() * g.cos()).pow(2) -
    FlatFac * h2 * (f.cos() * g.sin()).pow(2)
  )

proc rhoSinCosPhi*(geographLat, height: float64): tuple[rhoSinPhi, rhoCosPhi: float64] =
  ## Computes two quantities that are used elsewhere in the library
  let u = (geographLat.tan() * PolRad/EqRad).arctan()
  let x = height / (EqRad * 1000.0)

  let rhoSinPhi = (u.sin() * PolRad/EqRad) + (geographLat.sin() * x)
  let rhoCosPhi = u.cos() + (geographLat.cos() * x)

  (rhoSinPhi, rhoCosPhi)

proc rho*(geographLat: float64): float64 =
  ## Computes the distance from the Earth's center to a point on the
  ## ellipsoid
  0.9983271 +
    0.0016764 * (2.0 * geographLat).cos() -
    0.0000035 * (4.0 * geographLat).cos()

proc radOfParllLat*(geographLat: float64): float64 =
  ## Computes the radius of the parallel of a latitude
  EqRad * geographLat.cos() /
    (1.0 - (EccOfMeridian * geographLat.sin()).pow(2)).sqrt()

proc linearVelocityAtLat*(geographLat: float64): float64 {.inline.} =
  ## Computes the linear velocity of a point at a latitude
  RotAngularVelocity * radOfParllLat(geographLat)

proc radCurvOfMeridian*(lat: float64): float64 =
  ## Computes the radius of curvature of the Earth's meridian
  ## at a latitude
  EqRad * (1.0 - EccOfMeridian * EccOfMeridian) /
    (1.0 - (EccOfMeridian * lat.sin()).pow(2)).pow(1.5)

proc geographGeocentLatDiff*(geographLat: float64): float64 =
  ## Computes the difference between the geographic latitude and
  ## geocentric latitude
  angle.degFrmDMS(0, 0, 692.73) * (2.0 * geographLat).sin() -
    angle.degFrmDMS(0, 0, 1.16) * (4.0 * geographLat).sin()

proc equationOfTime*(
  JD: float64,
  sunAsc: float64,
  nutLong: float64,
  truOblq: float64
): float64 =
  ## Computes the equation of time in radians
  let t = time.julianMill(JD)
  let L = angle.limitTo360(
    280.4664567 +
    t * (360007.6982779 +
    t * (0.030328 +
    t * (1.0 / 49931.0 -
    t * (1.0 / 15300.0 +
    t / 2000000.0)))))
  
  (
    L - 0.0057183 - sunAsc.radToDeg() +
    nutLong.radToDeg() * truOblq.cos()
  ).degToRad()

proc anglBetwnDiurnalPathAndHz*(dec, observerLat: float64): float64 =
  ## Computes the angle between diurnal path and the horizon
  let B = dec.tan() * observerLat.tan()
  let C = (1.0 - B * B).sqrt()

  (C * dec.cos()).arctan2(observerLat.tan())
