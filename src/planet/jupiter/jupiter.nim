import angle
import coords
import math
import time
import planet
import nutation
import jupiter_moons
export jupiter_moons

proc eqSemidiameter*(jupEarthDist: float64): float64 {.inline.} =
  ## Computes Jupiter's equatorial semidiameter
  angle.degFrmDMS(0, 0, 98.44).degToRad() / jupEarthDist

proc polSemideaimeter*(jupEarthDist: float64): float64 {.inline.} =
  ## Computes Jupiter's polar semidiameter
  angle.degFrmDMS(0, 0, 92.06).degToRad() / jupEarthDist

type
  Ephemeris = object
    ##  Holds Jupiter's ephemeris values for physical observations
    De: float64
    Ds: float64
    P: float64
    w1: float64
    w2: float64

proc ephemeris*(
  JD: float64,
  mnOblq: float64,
  nutInLong: float64,
  nutInOblq: float64
): Ephemeris =
  ## Return quantites used in the ephemeris for physical observations
  ## of Jupiter
  let d = JD - 2433282.5
  let T1 = d / 36525.0

  let asc0 = (268.0 + 0.1061*T1).degToRad()
  let dec0 = (64.5  - 0.0164*T1).degToRad()

  let W1 = angle.limitTo360(17.710 + 877.90003539*d).degToRad()
  let W2 = angle.limitTo360(16.838 + 870.27003539*d).degToRad()

  let (l0, b0, R) = planet.heliocentCoords(planet.Planet.Earth, JD)

  var (l, b, r) = (0.0, 0.0, 0.0)
  var x = 0.0
  var y = 0.0
  var z = 0.0
  var jupEarthDist = 0.0
  var lightTime = 0.0

  var i: uint8 = 1
  while i <= 2:

    let (new_l, new_b, new_r) = planet.heliocentCoords(planet.Planet.Jupiter, JD - lightTime)
    l = new_l; b = new_b; r = new_r

    let (new_x, new_y, new_z) = planet.geocentEclRectCoords(l0, b0, R, l, b, r)
    x = new_x; y = new_y; z = new_z

    jupEarthDist = planet.distFrmEclRectCoords(x, y, z)
    lightTime = planet.lightTime(jupEarthDist)

    inc i

  l -= 0.01299.degToRad()*jupEarthDist / (r*r)
  (x, y, z) = planet.geocentEclRectCoords(l0, b0, R, l, b, r)
  jupEarthDist = planet.distFrmEclRectCoords(x, y, z)

  let asc_s = (mnOblq.cos()*l.sin() - mnOblq.sin()*b.tan()).arctan2(l.cos())
  let dec_s = (mnOblq.cos()*b.sin() + mnOblq.sin()*b.cos()*l.sin()).arcsin()

  let ds = (-dec0.sin()*dec_s.sin() - dec0.cos()*dec_s.cos()*(asc0 - asc_s).cos()).arcsin()

  let u = y*mnOblq.cos() - z*mnOblq.sin()
  let v = y*mnOblq.sin() + z*mnOblq.cos()
  var asc = u.arctan2(x)
  var dec = v.arctan2((x*x + u*u).sqrt())
  let zeta =
      (dec0.sin()*dec.cos()*(asc0 - asc).cos() - dec.sin()*dec0.cos())
      .arctan2(dec.cos()*(asc0 - asc).sin())

  let de = (-dec0.sin()*dec.sin() - dec0.cos()*dec.cos()*(asc0 - asc).cos()).arcsin()

  var w1 = angle.limitTo360(W1.radToDeg() - zeta.radToDeg() - 5.07033*jupEarthDist)
  var w2 = angle.limitTo360(W2.radToDeg() - zeta.radToDeg() - 5.02626*jupEarthDist)

  var C =
      57.2958 * (2.0*r*jupEarthDist + R*R - r*r - jupEarthDist*jupEarthDist) /
        (4.0 * r * jupEarthDist)
  if (l - l0).sin() < 0.0:
    C *= -1.0
  w1 = (w1 + C).degToRad()
  w2 = (w2 + C).degToRad()

  let truOblq = mnOblq + nutInOblq

  let q = 0.005693.degToRad()
  asc += q * (asc.cos()*l0.cos()*truOblq.cos() + asc.sin()*l0.sin()) / dec.cos()
  dec += q * (  l0.cos()*truOblq.cos()*(truOblq.tan()*dec.cos() -
                asc.sin()*asc.cos()) +
                asc.cos()*dec.sin()*l0.sin())

  let (ascNut, decNut) = nutation.nutationInEqCoords(
      coords.EqPoint(asc: asc, dec: dec),
      nutInLong,
      nutInOblq,
      truOblq
  )
  let asc1 = asc + ascNut
  let dec1 = dec + decNut

  let (asc0Nut, dec0Nut) = nutation.nutationInEqCoords(
      coords.EqPoint(asc: asc0, dec: dec0),
      nutInLong,
      nutInOblq,
      truOblq
  )
  let asc01 = asc0 + asc0Nut
  let dec01 = dec0 + dec0Nut

  let P = (dec01.cos() * (asc01 - asc1).sin())
          .arctan2(dec01.sin()*dec1.cos() - dec01.cos()*dec1.sin()*(asc01 - asc1).cos())

  Ephemeris(
    De: de,
    Ds: ds,
    P: P,
    w1: w1,
    w2: w2
  )
