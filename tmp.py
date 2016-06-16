#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import math
import numpy as np, matplotlib.pyplot as plt
import utils as ut
import dynamic as dy
from scipy.integrate import odeint
from pylab import *
from sympy import*


from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, h = symbols('f g h', cls=Function)


init_printing(use_unicode=True)
x, t, z, nu = symbols('x t z nu')

print(limit(sin(x)/x, x, 0))

print(latex(Integral(cos(x)**2, (x, 0, pi))))
# \int_{0}^{\pi} \cos^{2}{\left (x \right )}\, dx
plot(sin(2*sin(2*sin(x))))