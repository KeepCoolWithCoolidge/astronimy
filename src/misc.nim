import math

proc parllcAngl*(observerLat, hourAngle, dec: float64): float64 {.inline.} =
  ## Computes the parallactic angle of a celestial body
  hourAngle.sin().arctan2(
    observerLat.tan() * dec.cos() -
    hourAngle.cos() * dec.sin()
  )

proc parllcAnglOnHz*(observerLat, dec: float64): float64 {.inline.} =
  ## Computes the parallactic angle of a celestial body on the
  ## horizon
  (observerLat.sin() / dec.cos()).arccos()