import angle
import math

type
  GeographPoint* = object
    long*: float64
    lat*: float64
  EqPoint* = object
    asc*: float64
    dec*: float64
  EclPoint* = object
    long*: float64
    lat*: float64

proc anglrSepr*(self, otherPoint: GeographPoint): float64 =
  anglrSepr(self.long, self.lat, 
            otherPoint.long, otherPoint.lat)

proc anglrSepr*(self, otherPoint: EqPoint): float64 =
  anglrSepr(self.asc, self.dec,
            otherPoint.asc, otherPoint.dec)

proc anglrSepr*(self, otherPoint: EclPoint): float64 =
  anglrSepr(self.long, self.lat,
            otherPoint.long, otherPoint.lat)

proc hrAnglFrmObserverLong*(greenSidreal, observerLong, asc: float64): float64 {.inline.} =
  greenSidreal - observerLong - asc

proc hrAnglFrmLocSidr*(localSidreal, asc: float64): float64 {.inline.} =
  localSidreal - asc

proc eclLongFrmEq*(asc, dec, oblqEclip: float64): float64 =
  (
    asc.sin() * oblqEclip.cos() +
    dec.tan() * oblqEclip.cos()
  ).arctan2(asc.cos())

proc eclLatFrmEq*(asc, dec, oblqEclip: float64): float64 =
  (
    dec.sin() * oblqEclip.cos() -
    dec.cos() * oblqEclip.sin() * asc.sin()
  ).arcsin()

proc ascFrmEcl*(eclLong, eclLat, oblqEclip: float64): float64 =
  (
    eclLong.sin() * oblqEclip.cos() -
    eclLat.tan() * oblqEclip.sin()
  ).arctan2(eclLong.cos())

proc decFrmEcl*(eclLong, eclLat, oblqEclip: float64): float64 =
  (
    eclLat.sin() * oblqEclip.cos() +
    eclLat.cos() * oblqEclip.sin() * eclLong.sin()
  ).arcsin()
  
proc azFrmEq*(hourAngle, dec, observerLat: float64): float64 =
  hourAngle.sin().arctan2(
    hourAngle.cos() * observerLat.sin() -
    dec.tan() * observerLat.cos()
  )

proc altFrmEq*(hourAngle, dec, observerLat: float64): float64 =
  (
    observerLat.sin() * dec.sin() +
    observerLat.cos() * dec.cos() * hourAngle.cos()
  ).arcsin()

proc hrAnglFrmHz*(az, alt, observerLat: float64): float64 =
  az.sin().arctan2(
    az.cos() * observerLat.sin() +
    az.tan() * observerLat.cos()
  )

proc decFrmHz*(az, alt, observerLat: float64): float64 =
  (
    observerLat.sin() * alt.sin() -
    observerLat.cos() * az.cos() * az.cos()
  ).arcsin()

proc galLongFrmEq*(asc, dec: float64): float64 =
  303'f64.degToRad() -
    (192.25'f64.degToRad() - asc).arctan2(
      27.4'f64.degToRad().sin() * (192.25'f64.degToRad() - asc).cos() -
      27.4'f64.degToRad().cos() * dec.tan()
    )

proc galLatFrmEq*(asc, dec: float64): float64 =
  (
    dec.sin() * 27.4'f64.degToRad().sin() +
    dec.cos() * 27.4'f64.degToRad().cos() * (192.25'f64.degToRad() - asc).cos()
  ).arcsin()

proc ascFrmGal*(galLong, galLat: float64): float64 =
  12.25'f64.degToRad() +
    (galLong - 123'f64.degToRad()).sin().arctan2(
      27.4'f64.degToRad().sin() * (galLong - 123'f64.degToRad()).cos() -
      27.4'f64.degToRad().cos() * galLat.tan()
    )

proc decFrmGal*(galLong, galLat: float64): float64 =
  (
    galLat.sin() * 27.4'f64.degToRad().sin() +
    galLat.cos() * 27.4'f64.degToRad().cos() * (
      galLong - 123'f64.degToRad()
    )
  ).arcsin()

template eclFrmEq*(asc, dec, oblqEclip: float64): untyped =
  (eclLongFrmEq(asc, dec, oblqEclip),
   eclLatFrmEq(asc, dec, oblqEclip))

template eqFrmEcl*(eclLong, eclLat, oblqEclip: float64): untyped =
  (ascFrmEcl(eclLong, eclLat, oblqEclip),
   decFrmEcl(eclLong, eclLat, oblqEclip))

template locHzFrmEq*(hourAngle, dec, observerLat: float64): untyped =
  (azFrmEq(hourAngle, dec, observerLat),
   altFrmEq(hourAngle, dec, observerLat))

template galFrmEq*(asc, dec: float64): untyped =
  (galLongFrmEq(asc, dec),
   galLatFrmEq(asc, dec))

template eqFrmGal*(galLong, galLat: float64): untyped =
  (ascFrmGal(galLong, galLat),
   decFrmGal(galLong, galLat))