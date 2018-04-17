import angle
import math
import saturn_moons
export saturn_moons
import saturn_rings
export saturn_rings

const
  EquatorialUnitSemidiameter* = angle.degFrmDMS(0, 0, 82.73).degToRad()
  PolarUnitSemidiameter* = angle.degFrmDMS(0, 0, 73.82).degToRad()

proc apprntMagMuller*(delta, r, delU, B: float64): float64 =
  - 8.68 +
    5.0*(r*delta).log10() +
    0.044*delU.abs() -
    2.6*B.abs().sin() +
    1.25*B.sin().pow(2)

proc apprntMag84*(delta, r, delU, B: float64): float64 =
  - 8.88 +
    5.0*(r*delta).log10() +
    0.044*delU.abs() -
    2.6*B.abs().sin() +
    1.25*B.sin().pow(2)

proc polSemidiameter*(saturnEarthDist, earthLat: float64): float64 =
  let a = EquatorialUnitSemidiameter
  let b = PolarUnitSemidiameter
  let k = 1.0 - (b/a).pow(2)

  (a / saturnEarthDist) * (1.0 - k*earthLat.cos().pow(2)).sqrt()

proc eqSemidiameter*(saturnEarthDist: float64): float64 {.inline.} =
  EquatorialUnitSemidiameter / saturnEarthDist

