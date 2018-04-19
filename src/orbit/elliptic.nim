import angle
import math
from orbit import Node

proc trueAnom*(eccAnom, ecc: float64): float64 {.inline.} =
  2.0 * ((1.0 + ecc).sqrt() * (eccAnom / 2.0).tan()).arctan2(
    (1.0 - ecc).sqrt()
  )

proc radVecFrmEccAnom*(eccAnom, a, ecc: float64): float64 {.inline.} =
  a * (1.0 - ecc * eccAnom.cos())

proc radVecFrmTrueAnom*(trueAnom, a, ecc: float64): float64 {.inline.} =
  a * (1.0 - ecc * ecc) / (1.0 + ecc * trueAnom.cos())

proc eccAnom*(meanAnom, ecc, accuracy: float64): float64 {.inline.} =
  var prevE = 0.0
  var E = meanAnom

  while (E - prevE).abs() > accuracy:
    prevE = E
    E = meanAnom + ecc * E.sin()

  E

proc vel*(r, a: float64): float64 {.inline.} =
  42.1219 * (1.0 / r - 0.5 / a).sqrt()

proc perihVel*(a, e: float64): float64 {.inline.} =
  29.7847 * (
    (1.0 + e) / ((1.0 - e) * a)
  ).sqrt()

proc aphVel*(a, e: float64): float64 {.inline.} =
  29.7847 * (
    (1.0 - e) / ((1.0 + e) * a)
  ).sqrt()

proc lengthRamanujan*(a, b: float64): float64 {.inline.} =
  PI * (
    3.0 * (a + b) - ((a + 3.0 * b) * (3.0 * a + b)).sqrt()
  )

proc length*(a, b: float64): float64 {.inline.} =
  let A = (a + b) / 2.0
  let G = (a * b).sqrt()
  let H = (2.0 * a * b) / (a + b)

  PI * (21.0 * A - 2.0 * G - 3.0 * H) / 8.0

proc semimajAxis*(perih, ecc: float64): float64 {.inline.} =
  perih / (1.0 - ecc)

proc mnMotion*(semimajAx: float64): float64 {.inline.} =
  0.01720209895 / semimajAx.pow(1.5)

proc passThroughNode(
  v: float64,
  n: float64,
  a: float64,
  e: float64,
  T: float64
): (float64, float64) =
  let E = 2.0 * ((1.0 - e).sqrt() * (v / 2.0).tan()).arctan2((1.0 + e).sqrt())
  let M = E - e * E.sin()

  (T + M / n, a * (1.0 - e * E.cos()))

proc passageThroughNode*(
  w: float64,
  n: float64,
  a: float64,
  e: float64,
  T: float64,
  node: orbit.Node
): (float64, float64) =
  case node
  of orbit.Node.Ascend: passThroughNode(-w, n, a, e, T)
  of orbit.Node.Descend: passThroughNode(PI - w, n, a, e, T)

    