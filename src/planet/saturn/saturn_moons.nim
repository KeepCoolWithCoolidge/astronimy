import planet
import precess
import time
import math

type
  Moon* = enum
    Mimas
    Enceladus
    Tethys
    Dione
    Rhea
    Titan
    Hyperion
    Iapetus

type
  Info* = object
    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11: float64
    W0, W1, W2, W3, W4, W5, W6, W7, W8: float64
    s1, c1, s2, c2: float64
    e1: float64
    lambda0, beta0, delta: float64

proc createInfoStruct(JD: float64): Info =
  let angle1 = 28.0817.degToRad()
  let angle2 = 168.8112.degToRad()

  var info = Info(
    t1: JD - 2411093.0,
    t2: 0.0,
    t3: (JD - 2433282.423)/365.25 + 1950.0,
    t4: JD - 2411368.0,
    t5: 0.0,
    t6: JD - 2415020.0,
    t7: 0.0, t8: 0.0,
    t9: (JD - 2442000.5)/365.25,
    t10: JD - 2409786.0,
    t11: 0.0,
    W0: 0.0, W1: 0.0, W2: 0.0, W3: 0.0, W4: 0.0,
    W5: 0.0, W6: 0.0, W7: 0.0, W8: 0.0,
    s1: angle1.sin(),
    c1: angle1.cos(),
    s2: angle2.sin(),
    c2: angle2.cos(),
    e1: 0.0, lambda0: 0.0, beta0: 0.0, delta: 0.0
  )

  info.t2 = info.t1/365.25;
  info.t5 = info.t4/365.25;
  info.t7 = info.t6/36525.0;
  info.t8 = info.t6/365.25;
  info.t11 = info.t10/36525.0;

  info.W0 = 5.095 * (info.t3 - 1866.39).degToRad();
  info.W1 = (74.4     + 32.39*info.t2).degToRad();
  info.W2 = (134.3    + 92.62*info.t2).degToRad();
  info.W3 = (42.0     - 0.5118*info.t5).degToRad();
  info.W4 = (276.59   + 0.5118*info.t5).degToRad();
  info.W5 = (267.2635 + 1222.1136*info.t7).degToRad();
  info.W6 = (175.4762 + 1221.5515*info.t7).degToRad();
  info.W7 = (2.4891   + 0.002435*info.t7).degToRad();
  info.W8 = (113.35   - 0.2597*info.t7).degToRad();

  info.e1 = 0.05589 - 0.000346*info.t7;

  info

proc mimas(info: Info): (float64, float64, float64, float64) {.inline.} =
  let L = (
      127.64 +
      381.994497*info.t1 -
      43.57*info.W0.sin() -
      0.72*(3.0*info.W0).sin() -
      0.02144*(5.0*info.W0).sin()
  ).degToRad();

  let p = (106.1 + 365.549*info.t2).degToRad();
  let M = L - p;
  let C = (
      2.18287*M.sin() +
      0.025988*(2.0*M).sin() +
      0.00043*(3.0*M).sin()
  ).degToRad();

  let lambda1 = L + C;
  let gamma1 = 1.563.degToRad();
  let omega1 = (54.5 - 365.072*info.t2).degToRad();
  let r1 = 3.06879/(1.0 + 0.01905*(M + C).cos());

  (lambda1, gamma1, omega1, r1)

proc enceladus(info: Info): (float64, float64, float64, float64) {.inline.} =
  let L = (
      200.317 +
      262.7319002*info.t1 +
      0.25667*info.W1.sin() +
      0.20883*info.W2.sin()
  ).degToRad();

  let p = (309.107 + 123.44121*info.t2).degToRad();
  let M = L - p;
  let C = (
      0.55577*M.sin() +
      0.00168*(2.0*M).sin()
  ).degToRad();

  let lambda2 = L + C;
  let gamma2 = 0.0262.degToRad();
  let omega2 = (348.0 - 151.95*info.t2).degToRad();
  let r2 = 3.94118/(1.0 + 0.00485*(M + C).cos());

  result = (lambda2, gamma2, omega2, r2)

proc tethys(info: Info): (float64, float64, float64, float64) {.inline.} =
  let lambda3 = (
      285.306 +
      190.69791226*info.t1 +
      2.063*info.W0.sin() +
      0.03409*(3.0*info.W0).sin() +
      0.001015*(5.0*info.W0).sin()
  ).degToRad();
  let gamma3 = 1.0976.degToRad();
  let omega3 =(111.33 - 72.2441*info.t2).degToRad();
  let r3 = 4.880998;

  (lambda3, gamma3, omega3, r3)

proc dione(info: Info): (float64, float64, float64, float64) {.inline.} =
  let L = (
      254.712 +
      131.53493193*info.t1 -
      0.0215*info.W1.sin() -
      0.01733*info.W2.sin()
  ).degToRad();

  let p = (174.8 + 30.82*info.t2).degToRad();
  let M = L - p;
  let C = (
      0.24717*M.sin() +
      0.00033*(2.0*M).sin()
  ).degToRad();

  let lambda4 = L + C;
  let gamma4 = 0.0139.degToRad();
  let omega4 = (232.0 - 30.27*info.t2).degToRad();
  let r4 = 6.24871/(1.0 + 0.002157*(M + C).cos());

  (lambda4, gamma4, omega4, r4)

proc funroutine(e, a, omega, i, lambda1, p: float64; info: Info): (float64, float64, float64, float64) {.inline.} =
  let M = lambda1 - p;
  let C =   e*((2.0 - e*e*(0.25 - 0.0520833333*e*e))*M.sin() +
               e*((1.25 - 0.458333333*e*e)*(2.0*M).sin() +
                    e*((1.083333333 - 0.671875*e*e)*(3.0*M).sin() +
                         e*(1.072917*(4.0*M).sin() +
                             e*1.142708*(5.0*M).sin()))));
  let r = a*(1.0 - e*e)/(1.0 + e*(M + C).cos());
  let g = omega - 168.8112.degToRad();
  let a1 = i.sin()*g.sin();
  let a2 = info.c1*i.sin()*g.cos() - info.s1*i.cos();
  let gamma = (a1*a1 + a2*a2).sqrt().arcsin();
  let u = a1.arctan2(a2);
  let w = 168.8112.degToRad() + u;
  let h = info.c1*i.sin() - info.s1*i.cos()*g.cos();
  let phi = (info.s1*g.sin()).arctan2(h);
  let `lambda` = lambda1 + C + u - g - phi;

  (`lambda`, gamma, w, r)

proc rhea(info: Info): (float64, float64, float64, float64) {.inline.} =
  let p1 = (342.7 + 10.057*info.t2).degToRad();
  let a1 = 0.000265*p1.sin() + 0.01*info.W4.sin();
  let a2 = 0.000265*p1.cos() + 0.01*info.W4.cos();
  let e = (a1*a1 + a2*a2).sqrt();
  let p = a1.arctan2(a2);
  let N = (345.0 - 10.057*info.t2).degToRad();
  let lambda1 = (359.244 + 79.6900472*info.t1 + 0.086754*N.sin()).degToRad();
  let i = (28.0362 + 0.346898*N.cos() + 0.0193*info.W3.cos()).degToRad();
  let omega = (168.8034 + 0.736936*N.sin() + 0.041*info.W3.sin()).degToRad();
  let a = 8.725924;

  funroutine(e, a, omega, i, lambda1, p, info)

proc titan(info: Info): (float64, float64, float64, float64) {.inline.} =
  let L = (261.1582 + 22.57697855*info.t4 + 0.074025*info.W3.sin()).degToRad();
  let i1 = (27.45141 + 0.295999*info.W3.cos()).degToRad();
  let omega1 = (168.66925 + 0.628808*info.W3.sin()).degToRad();
  let a1 = info.W7.sin()*(omega1 - info.W8).sin();
  let a2 = info.W7.cos()*i1.sin() -
           info.W7.sin()*i1.cos()*(omega1 - info.W8).cos();
  let g0 = 102.8623.degToRad();
  let phi = a1.arctan2(a2);
  let s = (a1*a1 + a2*a2).sqrt();
  var g = info.W4 - omega1 - phi;
  var wDash = 0.0;

  var counter: uint8 = 1
  while counter <= 6: 
    wDash = info.W4 + 0.37515.degToRad()*((2.0*g).sin() - (2.0*g0).sin());
    g = wDash - omega1 - phi;
    inc counter

  let e1 = 0.029092 + 0.00019048*((2.0*g).cos() - (2.0*g0).cos());
  let q = 2.0*(info.W5 - wDash);
  let b1 = i1.sin()*(omega1 - info.W8).sin();
  let b2 = info.W7.cos()*i1.sin()*(omega1 - info.W8).cos() -
           info.W7.sin()*i1.cos();
  let theta = b1.arctan2(b2) + info.W8;
  let e = e1*(1.0 + 0.002778797*q.cos());
  let p = wDash + 0.159215.degToRad()*q.sin();
  let u = 2.0*(info.W5 - theta) + phi;
  let h = 0.9375*e1*e1*q.sin() + 0.1875*s*s*(2.0*(info.W5 - theta)).sin();
  let lambda1 = L - 0.254744.degToRad()*(info.e1*(info.W6.sin() + 0.75*info.e1*(2.0*info.W6).sin()) + h);
  let i = i1 + 0.031843.degToRad()*s*u.cos();
  let omega = omega1 + 0.031843.degToRad()*s*u.sin()/i1.sin();
  let a = 20.216193;

  funroutine(e, a, omega, i, lambda1, p, info)

proc hyperion(info: Info): (float64, float64, float64, float64) {.inline.} =
  let nu = (92.39 + 0.5621071*info.t6).degToRad();
  let et = (148.19 - 19.18*info.t8).degToRad();
  let theta = (184.8 - 35.41*info.t9).degToRad();
  let theta1 = theta - 7.5.degToRad();
  let a_ss = (176.0 + 12.22*info.t8).degToRad();
  let b_s = (8.0 + 24.44*info.t8).degToRad();
  let c_s = b_s + 5.0.degToRad();
  let wDash = (69.898 - 18.67088*info.t8).degToRad();
  let phi = 2.0*(wDash - info.W5);
  let xi = (94.9 - 2.292*info.t8).degToRad();
  let a = 24.50601 -
          0.08686*nu.cos() -
          0.00166*(et + nu).cos() +
          0.00175*(et - nu).cos();
  let e = 0.103458 -
          0.004099*nu.cos() -
          0.000167*(et + nu).cos() +
          0.000235*(et - nu).cos() +
          0.02303*et.cos() -
          0.00212*(2.0*et).cos() +
          0.000151*(3.0*et).cos() +
          0.00013*phi.cos();
  let p = wDash + (
      0.15648*xi.sin() -
      0.4457*nu.sin() -
      0.2657*(et + nu).sin() -
      0.3573*(et - nu).sin() -
      12.872*et.sin() +
      1.668*(2.0*et).sin() -
      0.2419*(3.0*et).sin() -
      0.07*phi.sin()
  ).degToRad();
  let lambda1 = (
      177.047 +
      16.91993829*info.t6 +
      0.15648*xi.sin() +
      9.142*nu.sin() +
      0.007*(2.0*nu).sin() -
      0.014*(3.0*nu).sin() +
      0.2275*(et + nu).sin() +
      0.2112*(et - nu).sin() -
      0.26*et.sin() -
      0.0098*(2.0*et).sin() -
      0.013*a_ss.sin() +
      0.017*b_s.sin() -
      0.0303*phi.sin()
  ).degToRad();
  let i = (
      27.3347 +
      0.643486*xi.cos() +
      0.315*info.W3.cos() +
      0.018*(theta.cos() - c_s.cos())
  ).degToRad();
  let omega = (
      168.6812 +
      1.40136*xi.cos() +
      0.68599*info.W3.sin() -
      0.0392*c_s.sin() +
      0.0366*theta1.sin()
  ).degToRad();

  funroutine(e, a, omega, i, lambda1, p, info)

proc iapetus(info: Info): (float64, float64, float64, float64) {.inline.} =
  let L = (261.1582 + 22.57697855*info.t4).degToRad();
  let w_dash1 = (91.796 + 0.562*info.t7).degToRad();
  let phi_big = (4.367 - 0.195*info.t7).degToRad();
  let theta = (146.819 - 3.198*info.t7).degToRad();
  let phi = (60.47 + 1.521*info.t7).degToRad();
  let pho = (205.055 - 2.091*info.t7).degToRad();
  let e1 = 0.028298 + 0.001156*info.t11;
  let w_dash0 = (352.91 + 11.71*info.t11).degToRad();
  let mu = (76.3852 + 4.53795125*info.t10).degToRad();
  let i1 = (
      18.4602 -
      info.t11*(0.9518 + info.t11*(0.072 - 0.0054*info.t11))
  ).degToRad();
  let omega1 = (
      143.198 -
      info.t11*(3.919 - info.t11*(0.116 + 0.008*info.t11))
  ).degToRad();
  let l = mu - w_dash0;
  let g = w_dash0 - omega1 - phi_big;
  let g1 = w_dash0 - omega1 - phi;
  let ls = info.W5 - w_dash1;
  let gs = w_dash1 - theta;
  let lT = L - info.W4;
  let gT = info.W4 - pho;
  let u1 = 2.0*(l + g - ls - gs);
  let u2 = l + g1 - lT - gT;
  let u3 = l + 2.0*(g - ls - gs);
  let u4 = lT + gT - g1;
  let u5 = 2.0*(ls + gs);
  let a = 58.935028 + 0.004638*u1.cos() + 0.058222*u2.cos();
  let e =
      e1 - 0.0014097*(g1 - gT).cos() +
      0.0003733*(u5 - 2.0*g).cos() +
      0.000118*u3.cos() + 0.0002408*l.cos() +
      0.0002849*(l + u2).cos() +
      0.000619*u4.cos();
  let w = (
      0.08077*(g1 - gT).sin() +
      0.02139*(u5 - 2.0*g).sin() -
      0.00676*u3.sin() +
      0.0138*l.sin() +
      0.01632*(l + u2).sin() +
      0.03547*u4.sin()
  ).degToRad();
  let p = w_dash0 + w/e1;
  let lambda1 = mu + (
      -0.04299*u2.sin() -
      0.00789*u1.sin() -
      0.06312*ls.sin() -
      0.00295*(2.0*ls).sin() -
      0.02231*u5.sin() +
      0.0065*(u5 + phi_big).sin()
  ).degToRad();
  let i = i1 + (
      0.04204*(u5 + phi_big).cos() +
      0.00235*(l + g1 + lT + gT + phi).cos() +
      0.0036*(u2 + phi).cos()
  ).degToRad();
  let w1 = (
      0.04204*(u5 + phi_big).sin() +
      0.00235*(l + g1 + lT + gT + phi).sin() +
      0.00358*(u2 + phi).sin()
  ).degToRad();
  let omega = omega1 + w1/i1.sin();

  funroutine(e, a, omega, i, lambda1, p, info)

proc D(X_j, Y_j, Z_j, D_j: float64; info: Info): tuple[X, Y, Z, D: float64] =
  let A1 = X_j;
  let B1 = info.c1*Y_j - info.s1*Z_j;
  let C1 = info.s1*Y_j + info.c1*Z_j;

  let A2 = info.c2*A1 - info.s2*B1;
  let B2 = info.s2*A1 + info.c2*B1;

  let A3 = A2*info.lambda0.sin() - B2*info.lambda0.cos();
  let B3 = A2*info.lambda0.cos() + B2*info.lambda0.sin();
  let C3 = C1;

  let A4 = A3;
  let B4 = B3*info.beta0.cos() + C3*info.beta0.sin();
  let C4 = C3*info.beta0.cos() - B3*info.beta0.sin();

  let et = A4;
  let nu = C4;

  let D = et.arctan2(nu);

  let X = A4*D_j.cos() - C4*D_j.sin();
  let Y = A4*D_j.sin() + C4*D_j.cos();
  let Z = B4;

  (X, Y, Z, D)

proc XYZ(lambdaJ, gammaJ, omegaJ, rJ: float64; info: Info, moon: Moon): tuple[X, Y, Z: float64] {.inline.} =
  let u = lambdaJ - omegaJ;
  let w = omegaJ - 168.8112.degToRad();

  ## moon of interest
  let X_j = rJ*(u.cos()*w.cos() - u.sin()*gammaJ.cos()*w.sin());
  let Y_j = rJ*(u.sin()*w.cos()*gammaJ.cos() + u.cos()*w.sin());
  let Z_j = rJ*u.sin()*gammaJ.sin();

  ## a ficticious ninth moon
  let X_9 = 0.0;
  let Y_9 = 0.0;
  let Z_9 = 1.0;

  ## some fancy stuff
  let (_, _, _, d9) = D(X_9, Y_9, Z_9, 0.0, info);
  var (X, Y, Z, _) = D(X_j, Y_j, Z_j, d9, info);

  ## correct for differential light-time
  let K = case moon
          of Mimas:     20947.0
          of Enceladus: 23715.0
          of Tethys:    26382.0
          of Dione:     29876.0
          of Rhea:      35313.0
          of Titan:     53800.0
          of Hyperion:  59222.0
          of Iapetus:   91820.0

  X += Z.abs()*(1.0 - (X/rJ).pow(2)).sqrt()/K;

  ## correct for the perspective effect
  let W = info.delta / (info.delta + Z/2475.0);
  X *= W;
  Y *= W;

  (X, Y, Z)


proc apprntRectCoords*(JD: float64, moon: Moon): tuple[X, Y, Z: float64] =
  var info = createInfoStruct(JD - 0.04942)

  let (planetEclPoint, saturnEarthDist) = planet.geocentApprntEclCoords(planet.Planet.Saturn, JD)
  var (lambda0, beta0) = (planetEclPoint.long, planetEclPoint.lat)

  (lambda0, beta0) = precess.precessEclCoords(
    lambda0, beta0,
    JD,
    time.julianDay(
      time.Date(
        year: 1950,
        month: time.Month.Jan,
        decimalDay: 1.5,
        calType: time.CalType.Gregorian
      )
    )
  )

  info.lambda0 = lambda0
  info.beta0 = beta0
  info.delta = saturnEarthDist

  let (lambdaJ, gammaJ, omegaJ, rJ) = case moon
                                      of Mimas: mimas(info)
                                      of Enceladus: enceladus(info)
                                      of Tethys: tethys(info)
                                      of Dione: dione(info)
                                      of Rhea: rhea(info)
                                      of Titan: titan(info)
                                      of Hyperion: hyperion(info)
                                      of Iapetus: iapetus(info)
  
  XYZ(lambdaJ, gammaJ, omegaJ, rJ, info, moon)