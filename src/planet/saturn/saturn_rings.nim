import angle
import coords
import planet
import time
import math

proc inc*(JC: float64): float64 {.inline.} =
  (
    28.075216 -
    JC * (
      0.012998 +
      JC * 0.000004
    )
  ).degToRad()

proc ascendNode*(JC: float64): float64 {.inline.} =
  (
    169.50847 +
    JC * (
      1.394681 +
      JC * 0.000412
    )
  ).degToRad()

type
  Elements = object
    B: float64
    B1: float64
    P: float64
    deltaU: float64
    a: float64
    b: float64

proc elements*(JD, nutInLong, trueOblq: float64): Elements =
  let (l0, b0, R) = planet.heliocentCoords(planet.Planet.Earth, JD)

  var (l, b, r) = (0.0, 0.0, 0.0)
  var (x, y, z) = (0.0, 0.0, 0.0)
  var saturnEarthDist = 0.0
  var lightTime = 0.0

  var i: uint8 = 1
  while i <= 2:
    let (new_l, new_b, new_r) =
      planet.heliocentCoords(planet.Planet.Saturn, JD - lightTime)
    l = new_l; b = new_b; r = new_r

    let (new_x, new_y, new_z) =
      planet.geocentEclRectCoords(l0, b0, R, l, b, r)
    x = new_x; y = new_y; z = new_z

    saturnEarthDist = planet.distFrmEclRectCoords(x, y, z)
    lightTime = planet.lightTime(saturnEarthDist)
    inc i
  
  let JC = time.julianCent(JD)
  let inc = inc(JC)
  let ascendNode = ascendNode(JC)

  var (`lambda`, beta) = planet.eclCoordsFrmEclRectCoords(x, y, z)
  let B = (
    inc.sin() * beta.cos() * (`lambda` - ascendNode).sin() -
    inc.cos() * beta.sin()
  ).arcsin()
  let semiMaj = angle.degFrmDMS(0, 0, 375.35).degToRad() / saturnEarthDist
  let semiMin = semiMaj * B.abs().sin()

  let N = (113.6655 + 0.8771 * JC).degToRad()

  let l1 = l - (0.01759 / r).degToRad()
  let b1 = b - (0.000764 * (l - N).cos() / r).degToRad()

  let B1 = (
    inc.sin() * b1.cos() * (l1 - ascendNode).sin() -
    inc.cos() * b1.sin()
  ).arcsin()

  let U1 = (
    inc.sin() * b1.sin() +
    inc.cos() * b1.cos() * (l1 - ascendNode).sin()
  ).arctan2(
    b1.cos() * (l1 - ascendNode).cos()
  )

  let U2 = (
    inc.sin() * beta.sin() +
    inc.cos() * beta.cos() * (`lambda` - ascendNode).sin()
  ).arctan2(
    beta.cos() * (`lambda` - ascendNode).cos()
  )
  let deltaU = (U1 - U2).abs()

  var lambda0 = ascendNode - 90.0.degToRad()
  let beta0 = 90.0.degToRad() - inc

  let q = 0.005693.degToRad()
  `lambda` += q * (l0 - `lambda`).cos() / beta.cos()
  beta += q * (l0 - `lambda`).sin() * beta.sin()

  lambda0 += nutInLong
  `lambda` += nutInLong

  let asc0 = coords.ascFrmEcl(lambda0, beta0, trueOblq)
  let dec0 = coords.decFrmEcl(lambda0, beta0, trueOblq)
  let asc = coords.ascFrmEcl(`lambda`, beta, trueOblq)
  let dec = coords.decFrmEcl(`lambda`, beta, trueOblq)

  let P = (dec0.cos() * (asc0 - asc).sin()). arctan2(
    dec0.sin() * dec.cos() - dec0.cos() * dec.sin() * (asc0 - asc).cos()
  )

  Elements(
    B: B,
    B1: B1,
    P: P,
    deltaU: deltaU,
    a: semiMaj,
    b: semiMin
  )

proc innEdgeOuterRing*(a, b: float64): (float64, float64) {.inline.} =
  (a * 0.8801, b * 0.8801)

proc outEdgeInnerRing*(a, b: float64): (float64, float64) {.inline.} =
  (a * 0.8599, b * 0.8599)

proc innEdgeInnerRing*(a, b: float64): (float64, float64) {.inline.} =
  (a * 0.665, b * 0.665)

proc innEdgeDuskRing*(a, b: float64): (float64, float64) {.inline.} =
  (a * 0.5486, b * 0.5486)
