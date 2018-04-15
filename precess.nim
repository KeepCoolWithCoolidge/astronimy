import angle
import coords
import time
import math

proc annualPrecess*(asc, dec, JD: float64): tuple[annPrecessAsc, annPrecessDec: float64] =
  ## Computes annual precession in equatorial coordinates towards a new
  ## epoch
  let JC = time.julianCent(JD)

  let m = (
    angle.degFrmHMS(0, 0, 3.07496) +
    angle.degFrmHMS(0, 0, 0.00186) * JC
  ).degToRad()

  let n = (
    angle.degFrmHMS(0, 0, 1.33621) -
    angle.degFrmHMS(0, 0, 0.00057) * JC
  ).degToRad()

  (m + n * asc.sin() * dec.tan(), n * asc.cos())

proc precessEqCoords*(
  oldAsc: float64,
  oldDec: float64,
  JD1: float64,
  JD2: float64
): tuple[newAsc, newDec: float64] =
  ## Computes equatorial coordinates reduced to a different epoch
  let T = time.julianCent(JD1)
  let t = (JD2 - JD1) / 36525.0

  let x = t * (angle.degFrmDMS(0, 0, 2306.2181) +
          T * (angle.degFrmDMS(0, 0, 1.39656) -
          T *  angle.degFrmDMS(0, 0, 0.000139)))
  let xi = (x + t * t * ((angle.degFrmDMS(0, 0, 0.30188) -
                      T * angle.degFrmDMS(0, 0, 0.000344)) +
                      t * angle.degFrmDMS(0, 0, 0.17998))).degToRad()
  let zeta = (x + t * t * ((angle.degFrmDMS(0, 0, 1.09468) -
                        T * angle.degFrmDMS(0, 0, 0.000066)) +
                        t * angle.degFrmDMS(0, 0, 0.018203))).degToRad()
  let y = T * angle.degFrmDMS(0, 0, 0.000217)
  let theta = (t * (angle.degFrmDMS(0, 0, 2004.3109) +
               T * (angle.degFrmDMS(0, 0, 0.8553) - y) -
              t * ((angle.degFrmDMS(0, 0, 0.42665) + y) +
                t * angle.degFrmDMS(0, 0, 0.041833)))).degToRad()
  let A = oldDec.cos() * (oldAsc + xi).sin()

  let B = theta.cos() * oldDec.cos() * (oldAsc + xi).cos() -
          theta.sin() * oldDec.sin()
  
  let C = theta.sin() * oldDec.cos() * (oldAsc + xi).cos() +
          theta.cos() * oldDec.sin()
  
  (A.arctan2(B) + zeta, C.arcsin())

proc precessEqCoordsFK5*(
  oldAsc: float64,
  oldDec: float64,
  JD1: float64,
  JD2: float64
): tuple[newAsc, newDec: float64] =
  ## Computes equatorial coordinates, from coordinates referred to the
  ## FK4 system, reduced to a different epoch
  let T = (JD1 - 2415020.3135) / 36524.2199
  let t = (JD2 - JD1) / 36524.2199

  let xi = t * (
    angle.degFrmDMS(0, 0, 2304.25) + T * angle.degFrmDMS(0, 0, 1.396) +
    t * (
      angle.degFrmDMS(0, 0, 0.302) + t * angle.degFrmDMS(0, 0, 0.018)
    )
  )

  let zeta = xi + t * t * (
    angle.degFrmDms(0, 0, 0.791) + t * angle.degFrmDMS(0, 0, 0.001)
  )

  let theta = t * (
    angle.degFrmDMS(0, 0, 2004.682) -
    T * angle.degFrmDMS(0, 0, 0.853) -
    t * (
      angle.degFrmDMS(0, 0, 0.426) +
      t * angle.degFrmDMS(0, 0, 0.042)
    )
  )

  let A = oldDec.cos() * (oldAsc + xi).sin()

  let B = theta.cos() * oldDec.cos() * (oldAsc + xi).cos() -
          theta.sin() * oldDec.sin()
  
  let C = theta.sin() * oldDec.cos() * (oldAsc + xi).cos() +
          theta.cos() + oldDec.sin()
  
  (A.arctan2(B) + zeta, C.arcsin())

proc anglesForEclChange(t, T: float64): tuple[nu, pi, rho: float64] {.inline.} = 
  let x = T * angle.degFrmDMS(0, 0, 0.000598)
  let nu = (
    angle.degFrmDMS(0, 0, 47.0029) -
    T * (angle.degFrmDMS(0, 0, 0.06603) - x) + 
    t * (
      (angle.degFrmDMS(0, 0, -0.03302) + x) +
      t * angle.degFrmDMS(0, 0, 0.00006)
    )
  ).degToRad()

  let pi = (
    174.876384 +
    T * (
      angle.degFrmDMS(0, 0, 3289.4789) +
      T * angle.degFrmDMS(0, 0, 0.60622)
    ) - 
    t * (
      (
        angle.degFrmDMS(0, 0, 869.8089) +
        T * angle.degFrmDMS(0, 0, 0.50491)
      ) -
      t * angle.degFrmDMS(0, 0, 0.03536)
    )
  ).degToRad()

  let y = T * angle.degFrmDMS(0, 0, 0.000042)
  let rho = (
    t * (
      angle.degFrmDMS(0, 0, 5029.0966) +
      T * (angle.degFrmDMS(0, 0, 2.22226) - y) +
      t * (
        (angle.degFrmDMS(0, 0, 1.11113) - y) -
        t * angle.degFrmDMS(0, 0, 0.000006)
      )
    )
  ).degToRad()

  (nu, pi, rho)

proc precessEclCoords*(
  oldLong: float64,
  oldLat: float64,
  jdOld: float64,
  jdNew: float64
): tuple[newLong, newLat: float64] =
  ## Computes ecliptic coordinates reduced to a different epoch
  let T = time.julianCent(jdOld)
  let t = (jdNew - jdOld) / 36525.0

  let (nu, pi, rho) = anglesForEclChange(t, T)

  let A = nu.cos() * oldLat.cos() * (pi - oldLong).sin() -
          nu.sin() * oldLat.sin()
  
  let B = oldLat.cos() * (pi - oldLong).cos()

  let C = nu.cos() * oldLat.sin() +
          nu.sin() * oldLat.cos() * (pi - oldLong).sin()
  
  let newLong = rho + pi - A.arctan2(B)
  let newLat = C.arcsin()

  (newLong, newLat)

proc precessOrbElements*(
  oldInc: float64,
  oldArgPerih: float64,
  oldLongAscendNode: float64,
  JD1: float64,
  JD2: float64
): tuple[newInc, newArgPerih, newLongAscendNode: float64] =
  ## Computes orbital elements reduced to a different equinox
  let T = time.julianCent(JD1)
  let t = (JD2 - JD1)/ 36525.0

  let (nu, pi, rho) = anglesForEclChange(t, T)

  var newInc: float64
  var newLongAscendNode: float64
  let phi = pi + rho

  if oldInc == 0.0:
    newInc = nu
    newLongAscendNode = phi + PI
  else:
    let A = oldInc.sin() * (oldLongAscendNode - pi).sin()

    let B = -nu.sin() * oldInc.cos() +
             nu.cos() * oldInc.sin() * (oldLongAscendNode - pi).cos()
    
    newInc = (A * A + B * B).sqrt().arcsin()
    newLongAscendNode = phi + A.arctan2(B)
  
  let deltaW = (-nu.sin() * (oldLongAscendNode - pi).sin() / newInc.sin()).arcsin()

  (newInc, oldArgPerih + deltaW, newLongAscendNode)