import angle
import coords
import math

proc eqHzParallax*(distToEarth: float64): float64 {.inline.} =
  (angle.degFrmDMS(0, 0, 8.794).degToRad().sin() / distToEarth).arcsin()

proc topocentEqCoords*(
  eqPoint: coords.EqPoint,
  eqHzParllx: float64,
  geographPoint: coords.GeographPoint,
  observerHt: float64,
  greenwSidr: float64
): coords.EqPoint
  let (rho_sin, rho_cos) = 