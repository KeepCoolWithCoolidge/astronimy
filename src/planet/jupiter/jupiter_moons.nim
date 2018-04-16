import angle
import coords
import math
import time

type
  Moon* = enum
    Io,
    Europa,
    Ganymede,
    Callisto

proc apprntRectCoords*(JD: float64, moon: Moon): tuple[X, Y: float64] =
  ## Computes the apparent rectangular coordinates for a Galilean moon
  let d = JD - 2451545.0
  let V = (172.74 + 0.00111588*d).degToRad()
  let M = (357.529 + 0.9856003*d).degToRad()
  let N = (20.02 + 0.0830853*d + 0.329*V.sin()).degToRad()
  let J = (66.115 + 0.9025179*d - 0.329*V.sin()).degToRad()
  let A = (1.915*M.sin() + 0.02*(2.0*M).sin()).degToRad()
  let B = (5.555*N.sin() + 0.168*(2.0*N).sin()).degToRad()
  let K = J + A - B
  let R = 1.00014 - 0.01671*M.cos() - 0.00014*(2.0*M).cos()
  let r = 5.20872 - 0.25208*N.cos() - 0.00611*(2.0*N).cos()
  let delta = (r*r + R*R - 2.0*r*R*K.cos()).sqrt()

  let phi = (R*K.sin()/delta).arcsin()

  let dMinusDeltaBy173 = d - delta/173.0
  let phiMinusB = phi - B

  let u1 = (163.8069 + 203.4058646*dMinusDeltaBy173).degToRad() + phiMinusB
  let u2 = (358.414  + 101.2916335*dMinusDeltaBy173).degToRad() + phiMinusB
  let u3 = (5.7176   + 50.234518*dMinusDeltaBy173).degToRad()   + phiMinusB

  var u = case moon
          of Io: u1
          of Europa: u2
          of Ganymede: u3
          of Callisto: (224.8092 + 21.48798*dMinusDeltaBy173).degToRad() + phiMinusB

  let G = (331.18 + 50.310482*dMinusDeltaBy173).degToRad()
  let H = (87.45  + 21.569231*dMinusDeltaBy173).degToRad()

  u = u + (case moon
           of Io: 0.473 * (2.0*(u1 - u2)).sin()
           of Europa: 1.065 * (2.0*(u2 - u3)).sin()
           of Ganymede: 0.165 * G.sin()
           of Callisto: 0.843 * H.sin()).degToRad()
  
  let rMoon = case moon
               of Io: (5.9057  - 0.0244*(2.0*(u1 - u2)).cos())
               of Europa: (9.3966  - 0.0882*(2.0*(u2 - u3)).cos())
               of Ganymede: (14.9883 - 0.0216*G.cos())
               of Callisto: (26.3627 - 0.1939*H.cos())

  let `lambda` = (34.35 + 0.083091*d + 0.329*V.sin()).degToRad() + B
  let ds = (3.12 * (`lambda` + 42.8.degToRad()).sin()).degToRad()
  let de = ds - (
      2.22*phi.sin()*(`lambda` + 22.0.degToRad()).cos() +
      1.3*(r - delta)*(`lambda` - 100.5.degToRad()).sin()/delta
  ).degToRad()

  let X =  rMoon * u.sin()
  let Y = -rMoon * u.cos() * de.sin()

  (X, Y)