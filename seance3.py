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
                    A4 = A[2:, 2:].copy()
                    B4 = B[2:].copy()
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


    return dico_vp, dico_vecteur_propre,dico_etat_lin4, dico_etat, dico_etat_lin



# *********************************************************************************************************************
# 5.2 - séance 3
# 5.2.1 - Vous choisirez un point de trim pour laquelle toutes les valeurs propres de la matrice A 4 sont
# à partie réelle strictement négative.



# On choisit le point de trim ci-dessous déjà utilisé lors de la séance 2
tout = v_propre(h_tuple, mac_tuple, km_tuple, ms_tuple)  # on récupere "tout"
vp = tout[0]  # on choisit uniquement les val propres, pas les vecteurs

# $$$ on fixe un pt de trim
pt_trim =  (3000.0, 0.5, 0.9, 1.)
vap1=vp[pt_trim] # on choisit les val propres d'un seul point de trim
# print(type(vap1))
vep = tout[1]  # on selectionne le dico des vecteurs propres
M = vep[pt_trim]  # on selectionne la matrice de vect propre associés au trim choisi ; M matrice de passage tq X = MX'
# X' vecteur d'état dans la base modale

A4B4_dict = tout[2]
A4B4 = A4B4_dict[pt_trim]
A4 = A4B4[0]
B4 = A4B4[1]
XU_dico = tout[3]
X,U = XU_dico[pt_trim]
X = np.array(X)
U=np.array(U)

# print(type(U))

#récup de A et B
AB_dict = tout[4]
AB = AB_dict[pt_trim]
A = AB[0]
B = AB[1]



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

#calcul du nouveau vecteur d'état intial on a alpha
#atan(wh/vae)
dalpha =  math.atan(2./X[dy.s_va])
Xi=X.copy()
Ui=np.copy(U)
Xi[dy.s_a]+=dalpha


print("*******************")

#X est l'état d'equilibre
dXi=np.zeros(6)
dXi[dy.s_a]=dalpha
A = np.array(A)
# def dyn_lin(dX,t):
#
#     dX = np.array(dX)
#     dXdot=np.dot(A,dX)
#     return dXdot

Xe = X.copy()



##########################################################################"
# point de trim différend de celui utilisé pour obtenir la matrice A et B
# X,U = XU_dico[(11000.0, 0.8, 0.9, 1.0)]
# X = np.array(X)
# Xe=X.copy()

# conclusion forte sensibilité aux conditions intiales

###########################################################################

def dyn_lin(X,t):
    dXdot=np.dot(A,X)-np.dot(A,Xe)
    return dXdot


def Xet(X,nb):
    '''
    produit
    :param X:
    :param nb:nb de ligne du tableau
    :return:un array contitué de nb lignes du vecteur X
    '''
    col = len(X)
    xet = np.zeros((nb,col))
    for i in range(nb):
        xet[i,:]= X
    return xet


def trajectoire(X, U, P):
    tt=240# 240s
    nb_val=1000
    t = np.linspace(0, tt, nb_val)
    sol = odeint(dy.dyn, X, t, (U, P))
    sol_lin = odeint(dyn_lin, X, t)

    dy.plot2(t, sol,sol_lin, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
    # dy.plot(t, soli, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
    plt.suptitle("Paramètres d'état en fonction du temps(secondes)- trim = {}".format(pt_trim))
    plt.savefig('seance3/param_of_time.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()
    plt.close()

    # plt.axis([0, max(sol[:, 0]), 10800, 11200])
    # plt.plot(sol[:, 0], sol[:, 1])
    # plt.grid(True)
    # plt.xlabel("Distance parcourue y (kilomètres) durant {}s".format(tt))
    # plt.ylabel("Altitude h (kilomètres)")
    # plt.title("Trajectoire d'un " + P.name)
    # plt.text(1, 10, "mouvement est rectiligne uniforme")
    # plt.text(1, 9, "Point de trim : mac = {:.2f}, h = {:.0f}m, ms = {:.1f}".format(mac, h, P.ms))
    # plt.text(1, 8, "masse = {:.0f}kg".format(P.m))
    # # plt.savefig('seance3/trajectoire.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    # plt.show()

# trajectoire(Xi,U,P)




# séance 3 :
# 5.2.4 - donner la repr modale

# modèle linéaire 4 composantes
dXi=np.zeros(4)
dXi[1]=dalpha

#calcul de M-1 l'inverse de M
M_1=np.linalg.inv(M)
dXip = np.dot(M_1,dXi)
print("*****************************************************************")
print("La matrice de passage est M = ")
print(affiche_mat(M))
print("*****************************************************************")
#vap1 contient les 4 val propres
diag4=np.diag(vap1) #np.exp(vap1)

print("Dans la base propre A est diagonale et vaut : ")
print(affiche_mat(diag4))
print()
print("La matrice de passage inverse M_1 = ")
print(affiche_mat(M_1) )
print()
B4p = np.dot(M_1,B4)
print("La matrice B4p = M_1*B4 dans la base propre est :")
print(affiche_mat(B4p))
print()
# on a A4 et B4 issus de la linéarisation autour de Xe

print("A4 = ", A4)
print()
print("B4 = ", B4)
# print("\n\n\n*****************************")

########################################################################
# Pour chaque pas de temps contenu dans t, on calcule les dXp
# tt=240
# t = np.linspace(0, tt, tt+1)
# ll=[]
# for time in t:
#     expDt = np.diag(np.exp(vap1*time))
#     dXp = np.dot(expDt,dXip)
#     dX=[0.,0.]
#     dX+=np.real(np.dot(M,dXp)).tolist()
#     Xt=X+np.array(dX)
#     # print(Xt)
#     ll.append(Xt)
# sol2=np.array(ll)
# # print(sol2)
# sol = odeint(dy.dyn, Xi, t, (U, P))
# dy.plot2(t, sol,sol2, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
# plt.show()

#dX = X(t)-Xe=X(t)-X => X(t) = X+dX

############################################################################################################"
# seconde méthode
# on calcule

# def dyn_lin(X, t, U, P):
#     '''
#     on implémente dXdot=AdX+BdU ; on a dU=0 Xdot=Xe + dXdot = Xdot + dXdot  X(t) = X+dX
#     return : Xdot
#     '''
#     dXdot = np.dot(A,dX)
#     Xdot= X+




