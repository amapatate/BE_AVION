#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import math
import numpy as np, matplotlib.pyplot as plt
import utils as ut
import dynamic as dy
from scipy.integrate import odeint
from pylab import *




markeur = ('o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', ' ')
plt.grid(True)
plt.axhline()
plt.axvline()
P = dy.Param_737_300()


h_tuple = (3000., 11000.)
h_kilometre = np.array(h_tuple) / 1000
mac_tuple = (0.5, 0.8)
km_tuple = (0.1, 0.9)
ms_tuple = (0.2, 1.)

ord_max = 1.
abs_min = 3
abs_max = 11
nb_ligne = len(km_tuple) * len(ms_tuple)


# on choisit le point de trim numéro 16 pour lequel on a :
h = 3000.
mac = 0.8
km = 0.9
ms = 1.

# intialisation de l'avion
P.set_mass_and_static_margin(km, ms)  # km masse avion, ms : marge statique

# calcul de va en m/s
va = dy.va_of_mach(mac, h)

# calcul de rho
rho = ut.isa(h)[1]

args = {}
args['h'] = h
args['va'] = dy.va_of_mach(mac, h)  # conv mach en m/s
X, U = dy.trim(P, args)
(A, B) = ut.num_jacobian(X, U, P, dy.dyn)


val_propre = np.linalg.eig(A)
print(np.real(val_propre[0]))
exit()

#calcul du nouveau vecteur d'état intial on a alpha
#atan(wh/vae)
dalpha =  math.atan(2./X[dy.s_va])
Xi=X.copy()
Xi[dy.s_a]+=dalpha


print("*******************")

#X est l'état d'equilibre
Xe = X.copy()

dXi=np.zeros(6)
dXi[dy.s_a]=dalpha


def dyn_lin(X,t):
    Xdot=np.dot(A,X)-np.dot(A,Xe)
    return Xdot

tt=240 # 240s
nb_val=tt+1
t = np.linspace(0, tt, nb_val)
sol_lin = odeint(dyn_lin, Xi, t)

dy.plot(t, sol_lin, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
plt.show()
print(sol_lin[240,:])

