import angle
import time
import math
import util

proc mnOblqLaskar*(JD: float64): float64 =
  ## Computes the mean obliquity of the ecliptic using
  ## J. Laskar's formula
  let u = time.julianCent(JD) / 100.0
  HornerEval(
    u,
    angle.degFrmDMS(23, 26, 21.448),
    [
      -angle.degFrmDMS(0, 0, 4680.93),
      -angle.degFrmDMS(0, 0, 1.55),
       angle.degFrmDMS(0, 0, 1999.25),
      -angle.degFrmDMS(0, 0, 51.38),
      -angle.degFrmDMS(0, 0, 249.67),
      -angle.degFrmDMS(0, 0, 39.05),
       angle.degFrmDMS(0, 0, 7.12),
       angle.degFrmDMS(0, 0, 27.87),
       angle.degFrmDMS(0, 0, 5.79),
       angle.degFrmDMS(0, 0, 2.45)
    ]
  ).degToRad()

proc mnOblqIAU*(JD: float64): float64 =
  ## Computes the mean obliquity of the ecliptic using
  ## the IAU formula
  let u = time.julianCent(JD) / 100.0
  HornerEval(
    u,
    angle.degFrmDMS(23, 26, 21.448),
    [
      -angle.degFrmDMS(0, 0, 46.815),
      -angle.degFrmDMS(0, 0, 0.00059),
       angle.degFrmDMS(0, 0, 0.001813)
    ]
  ).degToRad()

proc eclipPointsOnHz*(
  oblqEclip: float64, 
  observerLat: float64, 
  locSidreal: float64
): tuple[longPoint1, longPoint2: float64] =
  ## Computes the longitudes of the two ecliptic points on
  ## a horizon on Earth
  let p = (-locSidreal.cos()).arctan2(
    oblqEclip.sin() * observerLat.tan() +
    oblqEclip.cos() * locSidreal.sin()
  )

  (p, p + PI)

proc anglBetwnEclipAndHz*(
  oblqEclip: float64,
  observerLat: float64,
  locSidreal: float64
): float64 =
  ## Computes the angle between the ecliptic and a horizon
  ## on Earth
  (
    oblqEclip.cos() * observerLat.sin() -
    oblqEclip.sin() * observerLat.cos() * locSidreal.sin()
  ).arccos()

