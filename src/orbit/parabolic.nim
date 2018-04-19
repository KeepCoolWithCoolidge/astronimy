import angle
import math
from orbit import Node

proc trueAnomAndRadVec*(t, T, q: float64): (float64, float64) =
  let W = 0.03649116245 * (t - T) / q.pow(1.5)
  let G = W / 2.0
  let Y = (G + (G * G + 1.0).sqrt()).cbrt()
  let s = Y - 1.0 / Y
  let v = 2.0 * s.arctan()
  let r = q * (1.0 + s * s)

  (v, r)

proc passThroughNode(v, q, T: float64): (float64, float64) =
  let s = (v / 2.0).tan()
  let tNode = T + q.pow(1.5) * (s * (s * s + 3.0)) * 27.403895
  let radVec = q * (1.0 + s * s)

  (tNode, radVec)

proc passageThroughNode*(
  w: float64,
  q: float64,
  T: float64,
  node: orbit.Node
): (float64, float64) {.inline.} =
  case node
  of orbit.Node.Ascend: passThroughNode(-w, q, T)
  of orbit.Node.Descend: passThroughNode(PI - w, q, T)