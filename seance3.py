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
h = 11000.
mac = 0.8
km = 0.9
ms = 1.







# def etat_lin4(h_tuple, mac_tuple, km_tuple, ms_tuple):
#     '''
#     calcule les repr d'état linéarisées A et B
#     :param h_tuple: tuple contenant le altitudes
#     :param mac_tuple:
#     :param km_tuple:
#     :param ms_tuple:
#     :return:
#     '''
#     dico_etat = {}  # dico pour stocker les 16 états associés au 16 pts de trim
#     dico_etat_lin = {}
#     dico_etat_lin4 = {}
#
#     for i, km in enumerate(km_tuple):
#         for j, ms in enumerate(ms_tuple):
#             P.set_mass_and_static_margin(km, ms)
#             for idx, mac in enumerate(mac_tuple):
#                 dico = {}  # dico = args dans la fonction trim
#                 for h in h_tuple:
#                     dico['h'] = h
#                     dico['va'] = dy.va_of_mach(mac, h)  # conv mach en m/s
#                     # calcul des points de trim associés aux 4 valeurs h, mac, km et ms
#                     # et ajout dans le tuple etat
#                     etat = (X, U) = dy.trim(P, dico)
#                     etat_lin = (A, B) = ut.num_jacobian(X, U, P, dy.dyn)
#                     A4 = np.copy(A[2:, 2:])
#                     B4 = np.copy(B[2:])
#                     etat_lin4 = (A4, B4)
#                     # construction du tuple (h,mac,km,ms) qui jouera le role de clé du dico states
#                     pt_trim = (h, mac, km, ms)
#                     dico_etat[pt_trim] = etat
#                     dico_etat_lin[pt_trim] = etat_lin
#                     dico_etat_lin4[pt_trim] = etat_lin4
#     return dico_etat_lin4



# valeurs propres

def v_propre(h_tuple, mac_tuple, km_tuple, ms_tuple):
    '''

    :param h_tuple: tuple contenant le altitudes
    :param mac_tuple:
    :param km_tuple:
    :param ms_tuple:
    :return:
    '''
    dico_etat = {}  # dico pour stocker les 16 états associés au 16 pts de trim
    dico_etat_lin = {}
    dico_etat_lin4 = {}
    dico_vp = {}
    dico_vecteur_propre = {}

    for i, km in enumerate(km_tuple):
        for j, ms in enumerate(ms_tuple):
            P.set_mass_and_static_margin(km, ms)
            for idx, mac in enumerate(mac_tuple):
                dico = {}  # dico = args dans la fonction trim
                for h in h_tuple:
                    dico['h'] = h
                    dico['va'] = dy.va_of_mach(mac, h)  # conv mach en m/s
                    # calcul des points de trim associés aux 4 valeurs h, mac, km et ms
                    # et ajout dans le tuple etat
                    etat = (X, U) = dy.trim(P, dico)
                    etat_lin = (A, B) = ut.num_jacobian(X, U, P, dy.dyn)
                    A4 = np.copy(A[2:, 2:])
                    B4 = np.copy(B[2:])
                    etat_lin4 = (A4, B4)
                    # calcul vp val propre
                    vp = np.linalg.eig(A4)
                    # construction du tuple (h,mac,km,ms) qui jouera le role de clé du dico states
                    pt_trim = (h, mac, km, ms)
                    dico_etat[pt_trim] = etat
                    dico_etat_lin[pt_trim] = etat_lin
                    dico_etat_lin4[pt_trim] = etat_lin4
                    dico_vp[pt_trim] = vp[0]
                    dico_vecteur_propre[pt_trim] = vp[1:][0] # car vp est un tuple ; indice 1 => array des vect propre indice 0 val propres


    return dico_vp, dico_vecteur_propre,dico_etat_lin4, dico_etat



# *********************************************************************************************************************
# 5.2 - séance 3
# 5.2.1 - Vous choisirez un point de trim pour laquelle toutes les valeurs propres de la matrice A 4 sont
# à partie réelle strictement négative.



# On choisit le point de trim ci-dessous déjà utilisé lors de la séance 2
tout = v_propre(h_tuple, mac_tuple, km_tuple, ms_tuple)  # on récupere "tout"
vp = tout[0]  # on choisit uniquement les val propres, pas les vecteurs
pt_trim =  (11000.0, 0.8, 0.9, 1.0)
vap1=vp[pt_trim] # on choisit un seul point de trim

vep = tout[1]  # on selectionne le dico des vecteurs propres
M = vep[pt_trim]  # on selectionne la matrice de vect propre associés au trim choisi ; M matrice de passage tq X = MX'
# X' vecteur d'état dans la base modale

A4B4_dict = tout[2]
A4B4 = A4B4_dict[pt_trim]
A4 = A4B4[0]
B4 = A4B4[1]

toto = tout[3]
X,U = toto[pt_trim]
X = np.array(X)
U=np.array(U)
print(type(U))


# print(M)


def affiche_mat(M):
    try:
        nb_col=len(M[:,0])
        nb_ligne=len(M[0,:])
        for l in range(0, nb_ligne):
            for c in range(0, nb_col):
                f = "{:+.5f}  ".format(M[l,c])
                print(f,end='')
            print()
    except :
        for l in range(0, len(M)):
            f = "{:+.5f}  ".format(M[l])
            print(f)

#

# affiche_mat(M)
# print(A4)




# Pour la séance 3, le point de trim
# pt_trim =  (11000.0, 0.8, 0.9, 1.0)  a toutes ses vp à partie réelle <0
# valeurs propres  [
# -0.64246345 + 2.71414129j \
# -0.64246345 - 2.71414129j \
# -0.00269042 + 0.05730786j
# -0.00269042 - 0.05730786j
# ]

#5.2.1 - évolution sur 240s suite à une rafale de vent verticale de wh = 2 m/s induisant un nouvelle angle alpha
# initial à partir duquel on "lache" l'avion.

#calcul du nouveau vecteur d'état intial on a alpha = theta (incidence = assiette)
#atan(wh/vae)
dalpha =  math.atan(2./X[dy.s_va])
Xi=np.copy(X)
Ui=np.copy(U)
# Ui[dy.i_wz]+=2
Xi[dy.s_a]+=dalpha
Xi[dy.s_th]+=dalpha
print("*******************")
# affiche_mat(Xi)
# affiche_mat(X)
# print(dalpha)



def trajectoire(Xi,X, U, P):
    tt=240 # 240s
    t = np.linspace(0, tt, tt+1)
    sol = odeint(dy.dyn, X, t, (U, P))
    soli = odeint(dy.dyn, Xi, t, (U, P))

    dy.plot(t, sol, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
    dy.plot(t, soli, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
    plt.suptitle("Paramètres d'état en fonction du temps(secondes)- modèle non linéaire")
    plt.savefig('seance3/param_of_time.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()
    plt.close()

    plt.axis([0, max(sol[:, 0]) / 1000. + 10., 10.980, 11.020])
    plt.plot(sol[:, 0] / 1000., sol[:, 1] /1000.)
    plt.xlabel("Distance parcourue y (kilomètres) durant {}s".format(tt))
    plt.ylabel("Altitude h (kilomètres)")
    plt.title("Trajectoire d'un " + P.name)
    plt.text(1, 10, "mouvement est rectiligne uniforme")
    plt.text(1, 9, "Point de trim : mac = {:.2f}, h = {:.0f}m, ms = {:.1f}".format(mac, h, P.ms))
    plt.text(1, 8, "masse = {:.0f}kg".format(P.m))
    plt.savefig('seance3/trajectoire.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()

trajectoire(Xi,X,U,P)
