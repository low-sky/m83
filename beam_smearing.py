import numpy as np
import scipy.special as spec
import astropy.units as u
import astropy.constants as con

Md = 6.4e10 * u.M_sun
Rd = 2.9 * u.kpc

rad = np.linspace(0,10, 1001) * u.kpc

y = (rad / (2 * Rd)).to(u.dimensionless_unscaled)
y = y.value
V = (2 * con.G * Md / Rd * y**2 * (spec.iv(0, y) * spec.kv(0, y) -
                                   spec.iv(1, y) * spec.kv(1, y)))**(0.5)
V = V.to(u.km / u.s)
V = V.value
rad = rad.value

DvDr = np.diff(V[1:]) / (rad[1] - rad[0]) * 0.3 / 2
