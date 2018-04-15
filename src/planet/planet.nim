import angle
import coords
import time
import math
import VSOPD_87

type
  Planet* = enum
    Mercury
    Venus
    Earth
    Mars
    Jupiter
    Saturn
    Uranus
    Neptune

proc illumFracFrmPhaseAngl*(i: float64): float64 {.inline.} =
  (1.0 + i.cos()) / 2.0

proc illumFracFrmDist*(r, delta, R: float64): float64 {.inline.} =
  let x = r + delta
  (x * x - R * R) / (4.0 * r * delta)

proc phaseAngl*(r, delta, R: float64): float64 {.inline.} =
  (r * r + delta * delta - R * R) / (2.0 * r * delta)

proc posAngleOfBrightLimb*(
  sunEqPoint: coords.EqPoint,
  planetEqPoint: coords.EqPoint
): float64 =
  let a = sunEqPoint.dec.cos()
  let n = a * (sunEqPoint.asc - planetEqPoint.asc).sin()
  let d = sunEqPoint.dec.sin() * planetEqPoint.dec.cos() -
          planetEqPoint.dec.sin() * (sunEqPoint.asc - planetEqPoint.asc).cos() * a
  n.arctan2(d)

proc semidiameter*(
  planet: Planet,
  planetEarthDist: float64
): float64 =
  let s = case planet
          of Mercury: 
            angle.degFrmDMS(0, 0, 3.360).degToRad()
          of Venus: 
            angle.degFrmDMS(0, 0, 8.410).degToRad()
          of Mars:
            angle.degFrmDMS(0, 0, 4.680).degToRad()
          of Uranus:
            angle.degFrmDMS(0, 0, 35.02).degToRad()
          of Neptune:
            angle.degFrmDMS(0, 0, 33.50).degToRad()
          of Jupiter:
            #TODO: add jupiter type procedure
            0
          of Saturn:
            #TODO: add saturn type procedure
            0
          of Earth:
            #Add raise error
            0
  s / planetEarthDist

proc orbElements*(planet: Planet, JD: float64): tuple[L, a, e, i, omega, pi, M, w: float64] =
  let T = time.julianCent(JD)
  let TT = T * T
  let TTT = TT * T

  var L, a, e, i, omega, pi: float64

  case planet
  of Mercury:
    L = 252.250906 + 149474.0722491*T + 0.0003035*TT + 0.000000018*TTT
    a = 0.038709831
    e = 0.20563175 + 0.000020407*T - 0.0000000283*TT + 0.00000000018*TTT
    i = 7.004986 + 0.0018215*T - 0.0000181*TT + 0.000000056*TTT
    omega = 48.330893 + 1.1861883*T + 0.00017542*TT + 0.000000215*TTT
    pi = 77.456119 + 1.5564776*T + 0.00029544*TT + 0.000000009*TTT
  of Venus:
    L = 181.979801 + 58519.2130302*T + 0.00031014*TT + 0.000000015*TTT
    a = 0.72332982
    e = 0.00677192 - 0.000047765*T + 0.0000000981*TTT + 0.00000000046*TTT
    i = 3.394662 + 0.0010037*T - 0.00000088*TT - 0.000000007*TTT
    omega = 76.67992 + 0.9011206*T + 0.00040618*TT - 0.000000093*TTT
    pi = 131.563703 + 1.4022288*T - 0.00107618*TT - 0.000005678*TTT
  of Earth:
    L = 100.466457 + 36000.7698278*T + 0.00030322*TT + 0.00000002*TTT
    a = 1.000001018
    e = 0.01670863 - 0.000042037*T - 0.0000001267*TTT + 0.00000000014*TTT
    i = 0.0
    pi = 102.937348 + 1.7195366*T + 0.00045688*TT - 0.000000018*TTT
    omega = 0.0
  of Mars:
    L = 355.433 + 19141.6964471*T + 0.00031052*TT + 0.000000016*TTT
    a = 1.523679342
    e = 0.09340065 + 0.000090484*T - 0.0000000806*TTT - 0.00000000025*TTT
    i = 1.849726 - 0.0006011*T + 0.00001276*TT - 0.000000007*TTT
    omega = 49.558093 + 0.7720959*T + 0.00001557*TT - 0.000002267*TTT
    pi = 336.060234 + 1.8410449*T + 0.00013477*TT + 0.000000536*TTT
  of Jupiter:
    L = 34.351519 + 3036.3027748*T + 0.0002233*TT + 0.000000037*TTT
    a = 5.202603209 + 0.0000001913*T
    e = 0.04849793 + 0.000163225*T - 0.0000004714*TTT - 0.00000000201*TTT
    i = 1.303267 - 0.0054965*T + 0.00000466*TT - 0.000000002*TTT
    omega = 100.464407 + 1.0209774*T + 0.00040315*TT + 0.000000404*TTT
    pi = 14.331207 + 1.6126352*T + 0.00103042*TT - 0.000004464*TTT
  of Saturn:
    L = 50.077444 + 1223.5110686*T + 0.00051908*TT - 0.00000003*TTT
    a = 9.554909192 - 0.0000021390*T + 0.000000004*TT
    e = 0.05554814 - 0.000346641*T - 0.0000006436*TTT + 0.0000000034*TTT
    i = 2.488879 - 0.0037362*T - 0.00001519*TT + 0.000000087*TTT
    omega = 113.665503 + 0.877088*T - 0.00012176*TT - 0.000002249*TTT
    pi = 93.057237 + 1.9637613*T + 0.00083753*TT + 0.000004928*TTT
  of Uranus:
    L = 314.055005 + 429.8640561*T + 0.0003039*TT - 0.000000026*TTT
    a = 19.218446062 - 0.0000000372*T + 0.00000000098*TT
    e = 0.04638122 - 0.000027293*T + 0.0000000789*TTT + 0.00000000024*TTT
    i = 0.773197 + 0.0007744*T + 0.00003749*TT - 0.000000092*TTT
    omega = 74.005957 + 0.5211278*T + 0.00133947*TT + 0.000018484*TTT
    pi = 173.005291 + 1.486379*T + 0.00021406*TT + 0.000000434*TTT
  of Neptune:
    L = 304.348665 + 219.8833092*T + 0.00030882*TT + 0.000000018*TTT
    a = 30.110386869 - 0.0000001663*T + 0.00000000069*TT
    e = 0.00945575 + 0.000006033*T - 0.00000000005*TTT
    i = 1.769953 - 0.0093082*T - 0.00000708*TT + 0.000000027*TTT
    omega = 131.784057 + 1.1022039*T + 0.00025952*TT - 0.000000637*TTT
    pi = 48.120276 + 1.4262957*T + 0.00038434*TT + 0.00000002*TTT
  
  (
    angle.limitTo360(L).degToRad(),
    a,
    e,
    angle.limitTo360(i).degToRad(),
    angle.limitTo360(omega).degToRad(),
    angle.limitTo360(pi).degToRad(),
    angle.limitTo360(L - pi).degToRad(),
    angle.limitTo360(pi - omega).degToRad()
  )

proc heliocentCoords*(planet: Planet, JD: float64): tuple[long, lat, radVec: float64] =
  let VSOPD87_Terms = case planet
                      of Mercury: VSOPD_87.MercuryTerms
                      of Venus: VSOPD_87.VenusTerms
                      of Earth: VSOPD_87.EarthTerms
                      of Mars: VSOPD_87.MarsTerms
                      of Jupiter: VSOPD_87.JupiterTerms
                      of Saturn: VSOPD_87.SaturnTerms
                      of Uranus: VSOPD_87.UranusTerms
                      of Neptune: VSOPD_87.NeptuneTerms
  
  var L = 0.0
  var B = 0.0
  var R = 0.0

  let JM = time.julianMill(JD)

  var n: uint8 = 1
  for i in VSOPD87_Terms:
    var T = 1.0
    var y = 0.0
    for j in i:
      for k in j:
        y += k[0] * (k[1] + k[2] * JM).cos()
    case n
    of 1:
      L += y * T
    of 2:
      B += y * T
    of 3:
      R += y * T
    else:
      discard
    
    y = 0.0
    T *= JM

    inc n

  L = angle.limitToTwoPI(L)
  B = angle.limitToTwoPI(B)

  (L, B, R)

proc lightTime*(dist: float64): float64 {.inline.} =
  0.0057755183 * dist

proc geocentEclRectCoords*(L0, B0, R0, L, B, R: float64): tuple[X, Y, Z: float64] =
  let x = R * R.cos() * L.cos() - R0 * B0.cos() * L0.cos()
  let y = R * B.cos() * L.sin() - R0 * B0.cos() * L0.sin()
  let z = R * B.sin()           - R0 * B0.sin()

  (x, y, z)

proc eclCoordsFrmEclRectCoords*(x, y, z: float64): tuple[eclLong, eclLat: float64] {.inline.} =
  (
    y.arctan2(x),
    z.arctan2((x * x + y * y).sqrt())
  )

proc distFrmEclRectCoords*(x, y, z: float64): float64 {.inline.} =
  (x * y + y * y + z * z).sqrt()

proc geocentGeometEclCoords*(L0, B0, R0, L, B, R: float64): tuple[eclLong, eclLat, radVec, lightTime: float64] =
  let (x, y, z) = geocentEclRectCoords(L0, B0, R0, L, B, R)
  
  let (`lambda`, beta) = eclCoordsFrmEclRectCoords(x, y, z)
  let planetEarthDist = distFrmEclRectCoords(x, y, z)
  let lightTime = lightTime(planetEarthDist)

  (`lambda`, beta, planetEarthDist, lightTime)

proc geocentApprntEclCoords*(planet: Planet, JD: float64): tuple[planetEclPoint: coords.EclPoint, radVec: float64] =
  let (L0, B0, R0) = heliocentCoords(Earth, JD)

  let (L1, B1, R1) = heliocentCoords(planet, JD)
  let (_, _, _, t) = geocentGeometEclCoords(L0, B0, R0, L1, B1, R1)

  let (L2, B2, R2) = heliocentCoords(planet, JD - t)
  let (l2, b2, r2, _) = geocentGeometEclCoords(L0, B0, R0, L2, B2, R2)

  let eclPoint = coords.EclPoint(
    long: l2,
    lat: b2
  )

  (eclPoint, r2)

proc eclCoordsToFK5*(JD, eclLong, eclLat: float64): tuple[eclLongFK5, eclLatFK5: float64] =
  let JC = time.julianCent(JD)
  let lambda1 = eclLong - JC * (1.397 + JC * 0.00031).degToRad()
  let x = angle.degFrmDMS(0, 0, 0.03916).degToRad()

  let eclLongCorrection = - angle.degFrmDMS(0, 0, 0.09033).degToRad() +
                            x * (lambda1.cos() + lambda1.sin()) * eclLat.tan()
  (
    eclLong + eclLongCorrection,
    eclLat + x * (lambda1.cos() - lambda1.sin())
  )

proc geocentEqCoords*(
  X: float64,
  Y: float64,
  Z: float64,
  i: float64,
  w: float64,
  sigma: float64,
  oblqEclip: float64,
  v: float64,
  r: float64
): tuple[asc, dec, lightTime: float64] =
  let F = sigma.cos()
  let G = sigma.sin() * oblqEclip.cos()
  let H = sigma.sin() * oblqEclip.sin()

  let P = - sigma.sin() * i.cos()
  let Q =   sigma.cos() * i.cos() * oblqEclip.cos() - i.sin() * oblqEclip.sin()
  let R =   sigma.cos() * i.cos() * oblqEclip.sin() + i.sin() * oblqEclip.cos()

  let A = F.arctan2(P)
  let B = G.arctan2(Q)
  let C = H.arctan2(R)
  let a = (F*F + P*P).sqrt()
  let b = (G*G + Q*Q).sqrt()
  let c = (H*H + R*R).sqrt()

  let x = r * a * (A + w + v)
  let y = r * b * (B + w + v)
  let z = r * c * (C + w + v)

  let xi = X + x
  let nu = Y + y
  let et = Z + z

  let asc = angle.limitToTwoPI( nu.arctan2(xi) )
  let dec = et.arctan2( (xi*xi + nu*nu).sqrt() )
  let dist = (x*x + y*y + z*z).sqrt()

  (asc, dec, light_time(dist))

proc heliocentCoordsFrmOrbElements*(i, sigma, w, f, v, r: float64): tuple[long, lat: float64] =
  let u = w + v
  let x = r * (sigma.cos()*u.cos() - sigma.sin()*u.sin()*i.cos())
  let y = r * (sigma.sin()*u.cos() + sigma.cos()*u.sin()*i.cos())
  let z = r * i.sin() * u.sin()

  (y.arctan2(x), z.arctan2((x*x + y*y).sqrt()))

proc apprntMagMuller*(
  planet: Planet,
  i: float64,
  delta: float64,
  r: float64
): float64 =
  let x = 5.0 * (r * delta).log10()
  case planet
  of Mercury:
    x + 1.16 + (i - 50.0)*(0.02838 + (i - 50.0)*0.000102)
  of Venus:
    x - 4.0 + i*(0.01322 + i*i*0.0000004247)
  of Earth:
    raise newException(IOError, "Planet.Earth was passed to the procedure apprntMagMuller()")
  of Mars:
    x - 1.3 + i*0.01486
  of Jupiter:
    x - 8.93
  of Saturn:
    raise newException(IOError, "Planet.Saturn was passed to the procedure planet.apprntMagMuller()." &
      "Use the procedure planet.saturn.apprntMagMuller() instead.")
  of Uranus:
    x - 6.85
  of Neptune:
    x - 7.05

proc apprntMag84*(
  planet: Planet,
  i: float64,
  delta: float64,
  r: float64
): float64 =
  let x = 5.0 * (r * delta).log10()

  case planet
  of Mercury:
    x - 0.42 + i*(0.0380 - i*(0.000273 - i*0.00000200))
  of Venus:
    x - 4.40 + i*(0.0009 + i*(0.000239 - i*0.00000065))
  of Earth:
    #TODO: raise error
    0
  of Mars:
    x - 1.52 + i*0.016
  of Jupiter:
    x - 9.4 + i*0.005
  of Saturn:
    #TODO: raise error
    0
  of Uranus:
    x - 7.19
  of Neptune:
    x - 6.87

