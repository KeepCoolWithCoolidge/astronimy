import angle
import math

proc refracFrmApprntAlt15*(apprntAlt: float64): float64 =
  ## Computes the refraction term for true altitudes greater than 15
  ## degrees
  let x =
    angle.degFrmDMS(0, 0, 0.0668).degToRad() *
    (PI - apprntAlt).tan()
  
  angle.degFrmDMS(0, 0, 58.294).degToRad() * (PI - apprntAlt).tan() -
    x * x * x

proc refracFrmTrueAlt15*(trueAlt: float64): float64 =
  ## Computes the refraction term for apparent altitudes greater than 15
  ## degrees
  let x =
    angle.degFrmDMS(0, 0, 0.0824).degToRad() *
    (PI - trueAlt).tan()

  angle.degFrmDMS(0, 0, 58.276).degToRad() * (PI - trueAlt).tan() -
    x * x * x

proc refracFrmApprntAlt*(apprntAlt: float64): float64 =
  ## Computes the refraction term for apparent altitude
  result =  if apprntAlt == PI: 
              0.0
            else:
              let apprntAltDeg = apprntAlt.radToDeg()
              let a = apprntAltDeg + 7.31 / (apprntAltDeg + 4.4)
              let R = 1.0 / a.degToRad().tan()

              (R / 60.0).degToRad()

proc refracFrmTrueAlt*(trueAlt: float64): float64 =
  ## Computes the refraction term for true altitude
  result =  if trueAlt == PI:
              0.0
            else:
              let trueAltDeg = trueAlt.radToDeg()
              let a = trueAltDeg + 10.3 / (trueAltDeg + 5.11)
              let R = 1.02 / a.degToRad().tan()

              (R / 60.0).degToRad()

proc refracByPressr*(pressure: float64): float64 {.inline.} =
  ## Computes the refraction term modifier for local pressure (in millibars)
  pressure / 1010.0

proc refracByTemp*(temp: float64): float64 =
  ## Computes the refraction term modifier for local temperature (in Kelvin)
  283.0 / temp
