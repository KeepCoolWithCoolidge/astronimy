import angle
import orbit
import math

proc trueAnomAndRadVec*(
  t: float64,
  T: float64,
  ecc: float64,
  q: float64,
  accuracy: float64
): (float64, float64) =
  let daysFrmPerih = t - T
  if daysFrmPerih == 0.0:
    return (0.0, q)
  
  let k = consts.GAUSS_GRAV
  let d1 = 1000.0

  let q1 = k * ((1.0 + ecc) / q).sqrt() / (2.0 * q)
  let g = (1.0 - ecc) / (1.0 + ecc)

  let q2 = q1 * daysFrmPerih
  var s = 2.0 / (3.0 * q2.abs())
  s = 2.0 / (2.0 * (s.arctan() / 2.0).tan().cbrt().arctan()).tan()
  if daysFrmPerih < 0.0:
    s = -s
  if ecc != 1.0:
    var l = 0.0
    while true:
      let s0 = s
      var z = 1.0
      let y = s * s
      var g1 = -y * s
      var q3 = q2 * g * s * y * 2.0 / 3.0

      while true:
        z += 1.0
        g1 = -g1 * g * y
        let z1 = (z - (z + 1.0) * g) / (2.0 * z + 1.0)
        let f = z1 * g1
        q3 += f
        if z > 50.0 or f.abs() > d1:
          raise newException(IOError, "No convergence at near_parabolic.trueAnomAndRadVec")
        if f.abs() <= accuracy:
          break
      
      l += 1.0
      if l > 50.0:
        raise newException(IOError, "No convergence at near_parabolic.trueAnomAndRadVec")
      
      while true:
        let s1 = s
        s = (s * s * s * 2.0 / 3.0 + q3) / (s * s + 1.0)
        if (s - s1).abs() <= accuracy:
          break
      
      if (s - s0).abs() <= accuracy: break
    
  let v = angle.limitToTwoPi(2.0 * s.arctan())
  let r = q * (1.0 + ecc) / (1.0 + ecc * v.cos())

  (v, r)