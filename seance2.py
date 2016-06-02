__author__ = 'beldjiab'

# !/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np, matplotlib.pyplot as plt
import utils, dynamic

markeur = ('o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', ' ')
plt.grid(True)
P = dynamic.Param_737_300()

h_tuple = (3000., 11000.)
Ma_tuple = (0.5, 0.8)
ms_tuple = (0.2, 1.)
km_tuple = (0.1, 0.9)


P.set_mass_and_static_margin(0.1,0.2)
# dico = args dans la fonction trim

dico = {}
dico['h'] = 3000.
dico['va'] = dynamic.va_of_mach(0.5, 3000.) # mach et altitude h

print(dynamic.trim(P, dico))



