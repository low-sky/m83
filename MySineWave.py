# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt

def f(t):
    return np.sin(2*np.pi/l*t)
l=2.3   #l=wavelength
t=np.arange(-2*np.pi, 2*np.pi, 0.02)
plt.plot(t, f(t))
plt.ylabel('f(t)')
plt.xlabel('t')
plt.axis([-2*np.pi, 2*np.pi, -1, 1])
plt.show()
