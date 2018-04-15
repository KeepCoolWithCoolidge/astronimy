import angle
import math

proc mnAnnMotionOfCompan*(P: float64): float64 {.inline.} =
  ## Computes mean annual motion of companion star
  angle.TWO_PI / P

proc mnAnomOfCompan*(n, t, T: float64): float64 {.inline.} =
  ## Computes mean anomaly of companion star
  n * (t - T)

proc radVec*(a, e, eccAnom: float64): float64 {.inline.} =
  ## Computes radius vector of a binary star
  a * (1.0 - e * eccAnom.cos())

proc trueAnom*(e, eccAnom: float64): float64 {.inline.} =
  ## Computes true anomaly of a binary star
  2.0 * (
    ((1.0 + e) / (1.0 - e)).sqrt() * (eccAnom / 2.0).tan()
  ).arctan()

proc apprntCoordsAngl*(ascNodeCoords, trueAnom, w, i: float64): float64 =
  ## Computes apparent position angle of a binary sta
  let x = (
    (trueAnom + w).sin() * i.cos()
  ).arctan2((trueAnom + w).cos())

  angle.limitToTwoPi(x + ascNodeCoords)

proc anglrSepr*(radVec, trueAnom, w, i: float64): float64 =
  ## Computes angular separation of a binary star
  radVec * (
    ((trueAnom + w).sin() * i.cos()).pow(2) +
    (trueAnom + w).cos().pow(2)
  ).sqrt()

proc eccOfApprntOrb*(e, w, i: float64): float64 =
  ## Computes eccentricity of an apparent orbit
  let iCos = i.cos()
  let eWCos = e * w.cos()
  let eWCosSqr = eWCos * eWCos

  let a = (1.0 - eWCosSqr) * iCos * iCos
  let b = e * w.sin() * eWCos * iCos
  let c = 1.0 - eWCosSqr
  let d = ((a - c) * (a - c) + 4.0 * b * c).sqrt()

  ((2.0 * d) / (a + c + d)).sqrt()
  