import math

const
  TWO_PI*: float64 = 2.0 * PI

proc anglrSepr*(p1a1, p1a2, p2a1, p2a2: float64): float64 {.inline.}
proc degFrmDMS*(deg, min: int64; sec: float64): float64
proc dmsFrmDeg*(deg: float64): tuple[deg, min: int64; sec: float64] {.inline.}
proc hmsFrmDeg*(deg: float64): tuple[hour, min: int64; sec: float64] {.inline.}
proc limitTo360*(angle: float64): float64
proc limitToTwoPI*(angle: float64): float64

proc anglrSepr*(p1a1, p1a2, p2a1, p2a2: float64): float64 {.inline.} =
  ## Computes the angular separation between two angular points
  arccos(sin(p1a2) * sin(p2a2) + cos(p1a2) * cos(p2a2) * cos(p1a1 - p2a1))

proc degFrmDMS*(deg, min: int64; sec: float64): float64 =
  ## Computes an angle in degrees with decimals, from an angle
  ## expressed in degrees, arcminutes and arcseconds
  let (M, S) =
    if deg < 0: (-abs(min), -abs(sec))
    else: (min, sec)
  result = deg.float64 + M.float64 / 60.0'f64 + S / 3600.0'f64

proc dmsFrmDeg*(deg: float64): tuple[deg, min: int64; sec: float64] {.inline.} =
  ## Computes an angle expressed in degrees, arcminutes and
  ## arcseconds, from an angle in degrees with decimals
  result.deg = deg.int64
  let minutes = (deg - (result.deg.float64)) * 60.0'f64
  result.min = minutes.int64
  result.sec = (minutes - (result.min.float64)) * 60.0'f64

proc degFrmHMS*(hour, min: int64; sec: float64): float64 {.inline.} =
  15.0 * ((hour.float64) + (min.float64) / 60.0 + sec / 3600.0)

proc hmsFrmDeg*(deg: float64): tuple[hour, min: int64; sec: float64] {.inline.} =
  ## Computes an angle in degrees with decimals, from an angle
  ## expressed in hours, minutes and seconds
  let hours = deg / 15.0'f64
  result.hour = hours.int64
  let minutes = (hours - (result.hour.float64)) * 60.0'f64
  result.min = minutes.int64
  result.sec = (minutes - (result.min.float64)) * 60.0'f64

proc limitTo360*(angle: float64): float64 =
  ## Computes the equivalent angle in [0, 360] degree range
  let n = (angle / 360.0).int64
  let limitedAngle = angle - (TWO_PI * (n.float64))
  result = if limitedAngle < 0.0: limitedAngle + TWO_PI
           else: limitedAngle

proc limitToTwoPI*(angle: float64): float64 =
  ## Computes the equivalent angle in [0, 2π] radian range
  let n = (angle / TWO_PI).int64
  let limitedAngle = angle - (TWO_PI * (n.float64))
  result = if limitedAngle < 0.0: limitedAngle + TWO_PI
           else: limitedAngle
