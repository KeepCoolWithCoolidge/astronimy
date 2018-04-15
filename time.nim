import angle
import math

type
  CalType* = enum
    Gregorian
    Julian
  
  Month* = enum
    Jan = 1
    Feb
    Mar
    Apr
    May
    June
    July
    Aug
    Sept
    Oct
    Nov
    Dec

  Date* = object
    year*: int16
    month*: Month
    decimalDay*: float64
    calType*: CalType
  
  DayOfMonth* = object
    day*: range[1..31]
    hr*:  range[0..23]
    min*: range[0..59]
    sec*: range[0..60]
    timeZone*: float64
  
  Weekday* = enum
    Sunday
    Monday
    Tuesday
    Wednesday
    Thursday
    Friday
    Saturday

proc isLeapYear*(year: int16, calType: CalType): bool
proc julianDay*(date: Date): float64

proc weekdayFrmDate*(date: Date): Weekday =
  ## Computes `Weekday` from `Date`
  let date_OUT = Date(
    year: date.year, 
    month: date.month, 
    decimalDay: floor(date.decimalDay).float64,
    calType: Gregorian
    )
  let JD = julianDay(dateOut)
  let wd = (JD + 1.5).int64 mod 7
  
  case wd
  of 0:
    result = Sunday
  of 1:
    result = Monday
  of 2:
    result = Tuesday
  of 3:
    result = Wednesday
  of 4:
    result = Thursday
  of 5:
    result = Friday
  of 6:
    result = Saturday
  else:
    discard #Make this raise error

proc decimalDay*(day: DayOfMonth): float64 =
  ## Computes decimal day for a `DayOfMonth`
  result = (day.day.float64) +
           (day.hr.float64) / 24.0 +
           (day.min.float64) / 60.0 +
           day.sec.float64 / 60 -
           day.timeZone / 24.0

proc decimalYear*(date: Date): float64 =
  ## Computes decimal year for a `Date`
  var y = 0
  var days = 365.0
  var month = date.month.uint8

  if month > 1'u8: 
    inc(y, 31)
  if month > 2'u8:
    inc(y, 28)
    if isLeapYear(date.year, date.calType):
      inc(y)
      days += 1.0
  if month > 3'u8: 
    inc(y, 31)
  if month > 4'u8: 
    inc(y, 30)
  if month > 5'u8: 
    inc(y, 31)
  if month > 6'u8:
    inc(y, 30)
  if month > 7'u8:
    inc(y, 31)
  if month > 8'u8:
    inc(y, 31)
  if month > 9'u8:
    inc(y, 30)
  if month > 10'u8:
    inc(y, 31)
  if month > 11'u8:
    inc(y, 30)
  
  result = date.year.float64 + 
    (y.float64 + date.decimalDay) / days

proc isLeapYear*(year: int16, calType: CalType): bool =
  ## Checks if a year is a leap year
  result =  case calType
            of Julian:
              year mod 4 == 0
            of Gregorian:
              if year mod 100 == 0: 
                year mod 400 == 0
              else: 
                year mod 4 == 0

proc julianCent*(JD: float64): float64 {.inline.} =
  ## Computes Julian century for a Julian day 
  JD - 2451545.0 / 36525.0

proc julianMill*(JD: float64): float64 {.inline.} =
  ## Computes Julian millennium for a Julian day
  JD - 2451545.0 / 365250.0

proc julianDay*(date: Date): float64 =
  ## Computes Julian day from a `Date`
  let month = date.month.uint8
  let (y, m) =
    if month == 1'u8 or month == 2'u8:
      ((date.year - 1).float64, (month + 12).float64)
    else:
      (date.year.float64, month.float64)
  
  let a = floor(y / 100.0)
  let b = case date.calType
          of Gregorian:
            2.0 - a + floor(a / 4.0)
          of Julian:
            0.0
  
  result = floor(365.25 * (y + 4_716.0)) +
           floor(30.6001 * (m + 1.0)) +
           date.decimalDay +
           b -
           1_524.5

proc julianEphemerisDay*(JD, deltaT: float64): float64 {.inline.} =
  ## Computes the Julian Ephemeris day
  deltaT / 86_400.0 + JD

proc dateFrmJulianDay*(JD: var float64): tuple[year: int16, month: uint8, decimalDay: float64] =
  ## Computes a year, month and decimal day equivalent to a given Julian day
  assert(JD >= 0.0)
  JD += 0.5
  let Z = JD.int64
  let F = JD - Z.float64

  let A = if Z < 2_299_161:
            Z
          else:
            let alpha = floor(((Z.float64) - 1_867_216.25) / 36524.25).int64
            Z + 1 + alpha - floor(alpha.float64 / 4.0).int64
  
  let B = A + 1_524
  let C = floor((B.float64 - 122.1) / 365.25).int64
  let D = floor(365.25 * C.float64).int64
  let E = floor((B - D).float64 / 30.6001).int64
  
  let day = (B - D).float64 - floor(30.6001 * E.float64) + F
  let month = if E < 14:
                E - 1
              elif E == 14 or E == 15:
                E - 13
              else:
                0
  let year =  if month > 2:
                C - 4_716
              elif month == 1 or month == 2:
                C - 4_715
              else:
                0
  result = (year.int16, month.uint8, day)

proc apprntSidr*(mnSidr, nutInLong, trueOblq: float64): float64 =
  ## Computes apparent sidereal time from the mean sidereal time
  mnSidr + nutInLong * cos(trueOblq)

proc mnSidr*(JD: float64): float64 =
  ## Computes mean sidereal time for a Julian day
  let JC = julianCent(JD)
  angle.limitTo360(
    280.46061837 +
    360.98564736629 * (JD - 2_451_545.0) +
    JC * JC * (0.000387933 - JC / 387_100_000.0)
  ).radToDeg()

proc deltaT*(year: int32, month: uint8): float64 =
  ## Computes an approximate value of ΔT for a given year and month
  ## This function approximates ΔT from polynomial expressions using a  
  ## method different from that given in the *Meeus* book. The method
  ## used is given http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html;
  ## it covers a far wider time range, and is more accurate.
  let y = year.float64 + (month.float64 - 0.5) / 12.0

  if y < -500.0:
    let u = (y - 1_820.0) / 100.0
    return 32.0 * u * u - 20.0
  elif y < 500.0:
    let u = y / 100.0
    return 10_583.6 -
           u * (1_014.41 +
           u * (33.78311 +
           u * (5.952053 -
           u * (0.1798452 -
           u * (0.022174192 -
           u *  0.0090316521)))))
  elif y < 1_600.0:
    let u = (y - 1_000.0) / 100.0
    return 1_574.2 -
           u * (556.01 -
           u * (71.23472 +
           u * (0.319781 -
           u * (0.8503463 +
           u * (0.005050998 -
           u *  0.0083572073)))))
  elif y < 1_700.0:
    let u = y - 1_600.0
    return 120.0 -
           u * (0.9808 +
           u * (0.01532 -
           u /  7_129.0))
  elif y < 1_800.0:
    let u = y - 1_700.0
    return 8.83 +
           u * (0.1603 -
           u * (0.0059285 -
           u * (0.00013336 +
           u /  1_174_000.0)))
  elif y < 1_860.0:
    let u = y - 1_800.0
    return 13.72 +
           u * (0.332447 -
           u * (0.0068612 +
           u * (0.0041116 -
           u * (0.00037436 -
           u * (0.0000121272 +
           u * (0.0000001699 -
           u * (0.000000000875)))))))
  elif y < 1_900.0:
    let u = y - 1_860.0
    return 7.62 +
           u * (0.5737 -
           u * (0.251754 -
           u * (0.01680668 -
           u * (0.0004473624 +
           u /  233_174.0))))
  elif y < 1_920.0:
    let u = y - 1_900.0
    return -2.79 +
           u * (1.494119 -
           u * (0.0598939 -
           u * (0.0061966 +
           u *  0.000197)))
  elif y < 1_941.0:
    let u = y - 1_920.0
    return 21.20 +
           u * (0.84493 -
           u * (0.076100 -
           u * 0.0020936))
  elif y < 1_961.0:
    let u = y - 1_950.0
    return 29.70 +
           u * (0.407 -
           u * ((1.0 / 233.0) -
           u / 2_547.0))
  elif y < 1_986.0:
    let u = y - 1_975.0
    return 45.45 +
           u * (1.067 -
           u * ((1.0 / 260.0) +
           u / 718.0))
  elif y < 2_005.0:
    let u = y - 2_000.0
    return 63.86 +
           u * (0.3345 -
           u * (0.060374 -
           u * (0.0017275 +
           u * (0.000651814 +
           u * (0.0002373599)))))
  elif y < 2_050.0:
    let u = y - 2_000.0
    return 62.92 +
           u * (0.32217 +
           u *  0.005589)
  elif y <= 2_150.0:
    let u = (y - 1_820.0) / 100.0
    return 32.0 * u * u - 20.0 - 0.5628 * (2_150.0 - y)
  elif y > 2_150.0:
    let u = (y - 1_820.0) / 100.0
    return 32.0 * u * u - 20.0
  else:
    result = 0.0