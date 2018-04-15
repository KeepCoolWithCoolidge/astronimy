import math

proc brightnessRatio*(m1, m2: float64): float64 {.inline.} =
  ## Computes the brightness ratio of two stars
  10.pow(0.4 * (m2 - m1))

proc combinedMag*(m1, m2: float64): float64 {.inline.} =
  ## Computes the combined magnitude of two stars
  m2 - 2.5 * (brightnessRatio(m1, m2) + 1.0)

proc combinedMagOfMany*(m: openArray[float64]): float64 {.inline.} =
  ## Computes the combined magnitude of two or more stars
  var sum = 0.0
  for i in m:
    sum += 10.pow(-0.4 * i)
  -2.5 * sum.log10()

proc magDiff*(brightnessRatio: float64): float64 {.inline.} =
  ## Computes the difference in magnitude of two stars
  2.5 * brightnessRatio.log10()

proc absMagFrmParallax*(parallax, apparentMagnitude: float64): float64 {.inline.} =
  ## Computes the absolute magnitude of a star from its parallax
  apparentMagnitude + 5.0 + 5.0 * (parallax.radToDeg() * 3600.0).log10()

proc absMagFrmDist*(distance, apparentMagnitude: float64): float64 {.inline.} =
  ## Computes the absolute magnitude of a star from its distance from earth (parsecs)
  apparentMagnitude + 5.0 - 5.0 * distance.log10()

proc anglBetweenNorthCelesAndEclipPole*(eclipLong, eclipLat, oblqEclip: float64): float64 {.inline.} =
  ## Computes the angle between a vector from a star to the
  ## north celestial pole of the Earth and a vector from the
  ## same star to the north pole of the ecliptic
  (eclipLong.cos() * oblqEclip.tan()).arctan2(
    eclipLat.sin() * eclipLong.sin() * oblqEclip.tan() -
    eclipLat.cos()
  )

proc eqCoordsFrmMotion*(
  asc0: float64,
  dec0: float64,
  r: float64,
  deltaR: float64,
  properMotionAsc: float64,
  properMotionDec: float64,
  t: float64
): tuple[asc, dec: float64] =
  ## Computes the equatorial coordinates of a star at
  ## at a different time from it's motion in space
  ##
  ## This function Computes the equatorial coordinates
  ## of a star at a different time by taking into account
  ## it's proper motion, distance and radial velocity.
  let x = r * dec0.cos() * asc0.cos()
  let y = r * dec0.cos() * asc0.sin()
  let z = r * dec0.sin()

  let deltaAsc = 3600.0 * properMotionAsc.radToDeg() / 13751.0
  let deltaDec = 3600.0 * properMotionDec.radToDeg() / 206265.0

  let deltaX = (x / r) * deltaR - z * deltaDec * asc0.cos() - y * deltaAsc
  let deltaY = (y / r) * deltaR - z * deltaDec * asc0.sin() + x * deltaAsc
  let deltaZ = (z / r) * deltaR + r * deltaDec * dec0.cos()

  let x1 = x + t * deltaX
  let y1 = y + t * deltaY
  let z1 = z + t * deltaZ

  let asc = y1.arctan2(x1)
  let dec = z1.arctan2((x1 * x1 + y1 * y1).sqrt())

  (asc, dec)

proc properMotionInEqCoords*(
  asc: float64,
  dec: float64,
  pmotionAsc: float64,
  pmotionDec: float64,
  eclLat: float64,
  oblqEclip: float64
): tuple[pmotionLong, pmotionLat: float64] =
  let eclLatCos = eclLat.cos()

  let pmotionLong = (
    pmotionDec * oblqEclip.sin() * asc.cos() +
    pmotionAsc * dec.cos() * (
      oblqEclip.cos() * dec.cos() +
      oblqEclip.sin() * dec.sin() * asc.sin()
    )
  ) / (eclLatCos * eclLatCos)

  let pmotionLat = (
    pmotionDec * (
      oblqEclip.cos() * dec.cos() +
      oblqEclip.sin() * dec.sin() * asc.sin()
    ) -
    pmotionAsc * oblqEclip.sin() * asc.cos() * dec.cos()
  ) / eclLatCos

  (pmotionLong, pmotionLat)

