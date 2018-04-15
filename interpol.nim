import math
import util

proc threeValues*(y1, y2, y3, n: float64): float64 =
  ## Interpolates an intermediate value of a function from three of it's
  ## given values
  let a = y2 - y1
  let b = y3 - y2
  let c = b - a

  y2 + n * (a + b + n * c) / 2.0

proc fiveValues*(y1, y2, y3, y4, y5, n: float64): float64 =
  ## Interpolates an intermediate value of a function from five of it's
  ## given values
  let a = y2 - y1
  let b = y3 - y2
  let c = y4 - y3
  let d = y5 - y4

  let e = b - a
  let f = c - b
  let g = d - c

  let h = f - e
  let j = g - f

  let k = (j - h) / 12.0
  let h_j_12 = (h + j) / 6.0

  y3 +
    HornerEval(
      n,
      0.0,
      @[b + c - h_j_12,
      f - k,
      h_j_12,
      k]
    ) / 2.0
