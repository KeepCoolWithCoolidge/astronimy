import math

proc roundUptoDigits*(f: float64, decimalDigits: uint32): float64 =
  var d = 1.0
  
  for _ in 1..(decimalDigits + 1):
    d *= 10.0
  
  (f * d).round() / d

template HornerEval*(x, c: untyped, a: openarray[untyped]): untyped =
  var y = c
  var u = 1.0
  for item in a:
    u *= x
    y += u * item
  y