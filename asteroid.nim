import math

proc diameter*(absMag, albedo: float64): float64 {.inline.} =
  ## Computes the diameter (in meters) of an asteroid
  1_000.0 * 10'f64.pow(
    3.12 - absMag / 5.0 - 0.217147 * albedo.log10()
  )

proc apparentDiameter*(trueDiameter, asteroidEarthDist: float64): float64 {.inline.} =
  ## Computes the apparent diameter (in meters) of an asteroid
  1.3788 * trueDiameter / asteroidEarthDist
  