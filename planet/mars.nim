import angle
import coords
import math
import time
import planet

const
  NorthPolEqCoordsJ1950* = coords.EqPoint(
    asc: 317.342.degToRad(),
    dec: 52.7110.degToRad()
  )
  NorthPolEqCoordsJ2000* = coords.EqPoint(
    asc: 317.681.degToRad(),
    dec: 52.8860.degToRad()
  )

proc northPolEclCoords*(JC: float64): coords.EclPoint =
  coords.EclPoint(
    long: (352.9065 + 1.17330 * JC).degToRad(),
    lat:  (63.28180 - 0.00394 * JC).degToRad()
  )

type
  Ephemeris* = object
    De*: float64
    Ds*: float64
    P*: float64
    q*: float64
    w*: float64
    d*: float64

proc ephemeris*(
  JD: float64,
  northPoleEclCoords: coords.EclPoint,
  mnOblq: float64,
  nutInLong: float64,
  nutInOblq: float64
): Ephemeris =
  var (lambda0, beta0) = (northPoleEclCoords.long, northPoleEclCoords.lat)

  let (l0, b0, R) = planet.heliocentCoords(planet.Planet.Earth, JD)

  var (l, b, r) = (0.0, 0.0, 0.0)
  var (x, y, z) = (0.0, 0.0, 0.0)
  var marsEarthDist = 0.0
  var lightTime = 0.0
  
  var i: uint8 = 1
  while i <= 2:
    var (new_l, new_b, new_r) = planet.heliocentCoords(planet.Planet.Mars, JD - lightTime)
    l = new_l; b = new_b; r = new_r

    var (new_x, new_y, new_z) = planet.geocentEclRectCoords(l0, b0, R, l, b, r)
    x = new_x; y = new_y; z = new_z

    marsEarthDist = planet.distFrmEclRectCoords(x, y, z)
    lightTime = planet.lightTime(marsEarthDist)

    inc i
  
  var (`lambda`, beta) = planet.eclCoordsFrmEclRectCoords(x, y, z)

  let de = (
    -beta0.sin() * beta.sin() -
     beta0.cos() * beta.cos() * (lambda0 - `lambda`).cos()
  ).arcsin()

  let JC = time.julianCent(JD)
  let N = (49.5581 + 0.7721 * JC).degToRad()

  let l1 = l - (0.00697 / r).degToRad()
  let b1 = b - (0.000225 * (l - N).cos() / r).degToRad()
  let ds = (
    -beta0.sin() * b1.sin() -
     beta0.cos() * b1.cos() * (lambda0 - l1).cos()
  ).arcsin()

  let W = angle.limitTo360(
    11.504 +
    350.89200025 * (JD - lightTime - 2433282.5)
  ).degToRad()

  let asc0 = coords.ascFrmEcl(lambda0, beta0, mnOblq)
  let dec0 = coords.decFrmEcl(lambda0, beta0, mnOblq)

  let u = y*mn_oblq.cos() - z*mn_oblq.sin()
  let v = y*mn_oblq.sin() + z*mn_oblq.cos()
  let asc = u.arctan2(x)
  let dec = v.arctan2((x*x + u*u).sqrt())
  let zeta = (
      dec0.sin() * dec.cos() * (asc0 - asc).cos() -
      dec.sin() * dec0.cos()
  ).arctan2(dec.cos() * (asc0 - asc).sin())
  let w = W - zeta

  `lambda` += 0.005693.degToRad() * (l0 - `lambda`).cos()/beta.cos()
  beta += 0.005693.degToRad() * (l0 - `lambda`).sin() * beta.sin()

  `lambda` += nut_in_long
  lambda0 += nut_in_long
  let true_oblq = mn_oblq + nut_in_oblq

  let asc01 = coords.ascFrmEcl(lambda0, beta0, trueOblq)
  let dec01 = coords.decFrmEcl(lambda0, beta0, trueOblq)
  let asc1 = coords.ascFrmEcl(`lambda`, beta, trueOblq)
  let dec1 = coords.decFrmEcl(`lambda`, beta, trueOblq)

  let P = (dec01.cos() * (asc01 - asc1).sin()).arctan2(
      dec01.sin() * dec1.cos() - 
      dec01.cos() * dec1.sin() * (asc01 - asc1).cos()
  )

  let d =
      angle.degFrmDMS(0, 0, 9.36).degToRad() /
      marsEarthDist
  let k = planet.illumFracFrmDist(r, marsEarthDist, R)
  let q = (1.0 - k) * d

  Ephemeris(
    De: de,
    Ds: ds,
    P: P,
    q: q,
    w: w,
    d: d
  )