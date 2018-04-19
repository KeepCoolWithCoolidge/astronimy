const
  EQUATORIAL_RADIUS*: float64 = 6378137.0 ## Equatorial radius of the Earth (Semimajor axis) in meters
  FLATTENING*: float64 = 1.0 / 298.257223563 ## Flattening
  POLAR_RADIUS*: float64 = EQUATORIAL_RADIUS * (1.0 - FLATTENING) ## Polar radius of the Earth (Semiminor axis) in meters
  ANGULAR_VELOCITY*: float64 = 7.292115e-5 ## Angular velocity of the Earth in radians / second
  GRAV_CONST*: float64 = 3.986004418e+14 ## Gravitational constant (including the mass of the Earth's atmosphere) in meters^3 / second^2
