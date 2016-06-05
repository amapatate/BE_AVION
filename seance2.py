#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import math
import numpy as np, matplotlib.pyplot as plt
import utils as ut
import dynamic as dy
from scipy.integrate import odeint
from pylab import *

# import scipy as sc


markeur = ('o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', ' ')
plt.grid(True)
plt.axhline()
plt.axvline()
P = dy.Param_737_300()

# 4 Séance 2 : Linéarisation de la dynamique
# 4.1 Travail hors séances encadrées



h_tuple = (3000., 11000.)
h_kilometre = np.array(h_tuple) / 1000
mac_tuple = (0.5, 0.8)
km_tuple = (0.1, 0.9)
ms_tuple = (0.2, 1.)

ord_max = 1.
abs_min = 3
abs_max = 11
nb_ligne = len(km_tuple) * len(ms_tuple)


def trace(km, ms):
    P.set_mass_and_static_margin(km, ms)  # km masse avion, ms : marge statique
    for idx, mac in enumerate(mac_tuple):
        dico = {}  # dico = args dans la fonction trim
        alphas = []  # liste avec alpha1 obtenu avec h1=3000m et alpha2 associé à h2
        dphrs = []
        dths = []

        for h in h_tuple:
            dico['h'] = h
            dico['va'] = dy.va_of_mach(mac, h)  # conv mach en m/s
            X, U = dy.trim(P, dico)
            # alphas += [ut.deg_of_rad(X[3])]  # alpha en degrés
            # dphrs += [ut.deg_of_rad(U[0])] # delta PHR en degrés
            alphas += [X[3]]  # alpha en radians
            dphrs += [U[0]]  # delta PHR en radians
            dths += [U[1]]  # thrust de 0 à 1
        plt.title("ms = " + str(ms_tuple[0]) + "  km = " + str(km_tuple[0]) + " - " + P.name)
        # plt.title(r"$\alpha$")

        # grille 1x3 1 ligne 3 colonnes figure 1 : numérotation par ligne de gauche à droite
        plt.subplot(1, 3, 1)
        plt.plot(h_kilometre, alphas, marker=markeur[idx], label="Mac = " + str(mac))
        plt.plot(h_kilometre, alphas, marker=markeur[idx], label="Mac = " + str(mac))
        # plt.axis([abs_min, abs_max, 0., 20])
        plt.title(r"$\alpha$")

        # grille 1x3  figure 2
        plt.subplot(1, 3, 2)
        plt.subplots_adjust(wspace=0.4)
        plt.title(r"$\delta$phr")
        plt.plot(h_kilometre, dphrs, marker=markeur[idx], label="Mac = " + str(mac))
        # plt.axis([abs_min, abs_max, 0., 10])
        # plt.legend(loc=3,prop={'size':7})
        # plt.yticks(color='red')
        plt.xlabel('Altitude en km - ms = ' + str(ms) + " - coef. de masse km = " + str(km))  # , color='red')

        # grille 1x3  figure 3
        plt.subplot(1, 3, 3)
        plt.title(r"$\delta$th")
        plt.plot(h_kilometre, dths, marker=markeur[idx], label="Mac = " + str(mac))
        plt.legend(loc=1, prop={'size': 7})
        plt.axis([abs_min, abs_max, 0., ord_max])
        # plt.text(0.01 * abs_max, 0.8 * ord_max, "ms = " + str(ms_tuple[0]) + "  km = " + str(km_tuple[0]))
        # plt.legend(loc=2)
        # plt.xlabel('Altitude')
        # plt.ylabel(r'Incidence $\alpha$')


def trace_all(km_tuple, ms_tuple):
    for i, km in enumerate(km_tuple):
        for j, ms in enumerate(ms_tuple):
            trace(km, ms)
            plt.show()
            plt.close()


# trace_all(km_tuple,ms_tuple)


# 4.2.2 :
#
# 4.2.2.1 - Choisir arbitrairement un point de trim dans l’ensemble des points de trim proposés.

# 4.2.2.2 - Identifier l’équation de sustentation dans l’équation d’état.
# >>>>> c'est l'équation en "alpha point" on prend alpha point = q = 0 et alpha petit. on en déduit CL via la portance L



# 4.2.2.3 - En faisant l’hypothèse que l’incidence α est petite utiliser cette équation
# pour calculer la valeur du coefficient de la portance CL en fonction de la vitesse et de la masse.

# En tangage : Le compensateur le plus courant (quasiment indispensable) est le compensateur
# de profondeur, qui sert à équilibrer l'appareil sur l'axe de tangage (tendance à piquer ou à cabrer).
# Pour une masse donnée, à chaque vitesse de vol correspond une incidence de vol et donc une position
#  de la gouverne de profondeur. À chaque changement de vitesse, l'incidence change,
# ainsi il faut compenser à nouveau (ou trimer) pour garder cette nouvelle incidence.
# Le trim correspond à une incidence et donc une vitesse, pour une masse donnée de l'avion.




# on choisit le point de trim numéro 16 pour lequel on a :
h = 11000.
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
alpha16 = X[3]  # alpha en radians
dphr16 = U[0]  # delta PHR en radians
dth16 = U[1]  # thrust de 0 à 1

# calcul de la poussée pour le point de trim 16
F16 = dy.propulsion_model(X, U, P)

# calcul de CL16, la val du CL au point de trim16
CL16 = 2 * (P.m * P.g - math.sin(alpha16)) / ((rho * np.power(va, 2) * P.S))

print("valeurs obtenues via la fonction trim.")
print("alpha16 = ", alpha16, " rad")
print("angle réglage PHR dphr16 = ", dphr16, " rad")
print("Gaz ou thrust dth16 = ", dth16)
print("CL16 = ", CL16)

# 4.2.2.4 Conclure quant à la technique de réglage de la vitesse de l’avion.

# >>>> pour augmenter la vitesse, il faut diminuer CL donc alphae donc augmenter algébriquement
# l'angle du PHR pour augmenter la portance ou diminuer la déportance.

# 4.2.2.5 - Utiliser la valeur de CL obtenue ainsi que les graphiques tracés lors de la première séance
# pour déterminer (approximativement) les valeurs de trim αe et δPHRe .

print("La valeur de alphae obtenue via le graphique de la séance 1 CLe=f(alphae) est ", ut.rad_of_deg(6.7))
print(
    "La valeur de dphre obtenue via le graphique dphre = f(alphae) de la séance 1 en utilisant la valeur de alphae précédente est -0.21rad")

# 4.2.2.6 - Comment obtenir la valeur de trim de la manette des gaz δth à partir de la polaire équilibrée ?
# connaissant CL16=CLe = 0.57 la polaire équilibrée donne CDe = 0.038 donc f = CLe/CDe donc D = mg/f ;
# on utilise F = Fmax * delta-th = D, i.e. poussee = trainee, pour calculer deltath
# le graphe de poussée donne Fmax = 50128N
# delta_th = D/Fmax
# par le calcul on delta_th = 0.67 ; par la méthode graphique on a delta_th=0.0.70


CLe = CL16
CDe = 0.038
f = CLe / CDe
# calcul du coef de poussée tel que F = coef * deltath
coef = dy.coef_pousseel(X, U, P)

# autre méthode
Fmax = 50128.

# calcul de la trainée
D = dy.trainee(h, mac, CDe, P)

# deltath
deltath_polaire = D / Fmax

print("deltath obtenu via la polaire équilibrée = ", deltath_polaire)

# autre méthode de calcul de delta_th : D = mg/f  ; f = L/D=CLe/CDe
print("D1 = ", D)
D = P.g * P.m / f
print("D2 = ", D)
print("autre méthode f = CLe/CDe delta_th = ", D / Fmax)

# 4.2.3 :
# 4.2.3.1 :
# Pour le point de trim précédemment étudié vérifier que le calcul numérique conduit aux
# mêmes résultats que la méthode graphique.
# 4.2.3.2 :
# simuler la trajectoire de l’avion sur 100s.
# le vol est en palier avec les points de trim alpha_e, delta_phre et delta_the
# on a q = 0 d'où theta = cte. Or, en palier gamma=0, donc theta = alpha_e = cte
# va = cte q_point = 0 M =0
# D = F * cos(alpha)
# L + F * sin(alpha) = mg
# il reste l'ODE y_point = va  d'où y = va * t

# dy.plot(np.linspace(0,100,20),np.array(X))

print("******************************************************************")
print('')
# print(X)
print("******************************************************************")
print('')


# print(U)



def trajectoire(X, U, P):
    t = np.linspace(0, 100, 101)
    sol = odeint(dy.dyn, X, t, (U, P))
    dy.plot(sol[:, 0], sol, U=None, figure=None, window_title="Paramètres d'état en fonction de y(mètres)")
    plt.savefig('seance2/param_of_y.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()
    plt.close()

    dy.plot(t, sol, U=None, figure=None, window_title="Paramètres d'état en fonction du temps(secondes)")
    plt.savefig('seance2/param_of_time.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()
    plt.close()

    plt.axis([0, max(sol[:, 0]) / 1000. + 10., 0, 12])
    plt.plot(sol[:, 0] / 1000., sol[:, 1] / 1000.)
    plt.xlabel("Distance parcourue y (kilomètres) durant 100s")
    plt.ylabel("Altitude h (kilomètres)")
    plt.title("Trajectoire d'un " + P.name)
    plt.text(1, 10, "Comme attendu le mouvenement est rectiligne uniforme")
    plt.text(1, 9, "Point de trim : mac = {:.2f}, h = {:.0f}m, ms = {:.1f}".format(mac, h, P.ms))
    plt.text(1, 8, "masse = {:.0f}kg".format(P.m))
    plt.savefig('seance2/trajectoire.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()


# trajectoire(X, U, P)


# 4.2.4

# 4.2.4.1 : Linéariser numériquement le modèle d’état pour toutes les conditions de trim.
#


# 4.2.4.2 Extraire des représentations d’état linéarisées {A,B} de dimension 6
# (qui est la dimension du vecteur d'état y h  va α θ q),les matrices A4 et B4 associées aux composantes va α θ q
# du vecteur d’état.
# A4=A[2:,2:] B4 = B[2:]


# 4.2.4.3 Calculer numériquement les valeurs propres de la matrice A4 .
#
# 4.2.4.4 : Vous tracerez les valeurs propres pour toutes les conditions de trim en
# faisant varier successivement l’altitude h, le nombre de Mach Ma , la marge statique ms
# et la masse m (au travers les différentes valeurs du coefficient km de réglage de la masse).

# 4.2.4.5 Comment varient les valeurs propres en fonction de ces différents paramètres ?


# 4.2.4.1 :

def etat_lin4(h_tuple, mac_tuple, km_tuple, ms_tuple):
    '''
    calcule les repr d'état linéarisées A et B
    :param h_tuple: tuple contenant le altitudes
    :param mac_tuple:
    :param km_tuple:
    :param ms_tuple:
    :return:
    '''
    dico_etat = {}  # dico pour stocker les 16 états associés au 16 pts de trim
    dico_etat_lin = {}
    dico_etat_lin4 = {}

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
                    A4 = A[2:, 2:]
                    B4 = B[2:]
                    etat_lin4 = (A4, B4)
                    # construction du tuple (h,mac,km,ms) qui jouera le role de clé du dico states
                    pt_trim = (h, mac, km, ms)
                    dico_etat[pt_trim] = etat
                    dico_etat_lin[pt_trim] = etat_lin
                    dico_etat_lin4[pt_trim] = etat_lin4
    return dico_etat_lin4


# toto = etat_lin4(h_tuple, mac_tuple, km_tuple, ms_tuple)
# h, mac, km, ms = 11000, 0.8, 0.9, 1.
# print(toto[(h, mac, km, ms)][1])  # indice 0 donne A, l'incdice 1 donne B pour le pt de trim spécifié par le quadruplet


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
                    A4 = A[2:, 2:]
                    B4 = B[2:]
                    etat_lin4 = (A4, B4)
                    # calcul vp val propre
                    vp = np.linalg.eig(A4)
                    # construction du tuple (h,mac,km,ms) qui jouera le role de clé du dico states
                    pt_trim = (h, mac, km, ms)
                    dico_etat[pt_trim] = etat
                    dico_etat_lin[pt_trim] = etat_lin
                    dico_etat_lin4[pt_trim] = etat_lin4
                    dico_vp[pt_trim] = vp[0]
                    dico_vecteur_propre[pt_trim] = vp[1:]

    return dico_vp, dico_vecteur_propre


vp = v_propre(h_tuple, mac_tuple, km_tuple, ms_tuple)[0]  # on choisit uniquement les val propres, pas les vecteurs


def affiche_vp(h_tuple, mac_tuple, km_tuple, ms_tuple,P):
    num_plot = 1
    fig = plt.figure(1)
    for idx,h in enumerate(h_tuple):
        if num_plot > 8:
            num_plot = 1
            plt.show()
            plt.close()
        for km in km_tuple:
            for ms in ms_tuple:
                for mac in mac_tuple:
                        pt_trim = (h, mac, km, ms)
                        vpp = vp[pt_trim]
                        x = [vpp[k2].real for k2 in range(0, len(vpp))]
                        y = [vpp[k2].imag for k2 in range(0, len(vpp))]
                        # for i,xx in enumerate(x):
                        #         plt.annotate("toto",xy=(xx, y[i]), xytext=(xx-0.2, y[i]-0.2),
                        #                  arrowprops={'facecolor': 'red', 'shrink': 0.05})
                        plt.subplot(4, 2, num_plot)
                        plt.scatter(x, y, color='red',s=4)
                        # plt.xscale('symlog')
                        # plt.yscale('symlog')
                        plt.grid(True)
                        plt.axis([-3,0.1,-6.3,6.3])
                        plt.text(-2.5,4,"Ma={:.1f}, km = {:.1f} , ms = {:.1f}".format(pt_trim[1],pt_trim[2],pt_trim[3]), fontsize=8)
                        fig.subplots_adjust(hspace=0.3)
                        plt.suptitle("Valeurs propres pour h = {:.0f} km pour le {}".format(h_tuple[idx]/1000,P.name))
                        print("pt_trim = ", pt_trim, " \n valeurs propres ", vpp, "\n")
                        plt.xticks(fontsize=7)
                        plt.yticks(fontsize=7)
                        plt.axhline()
                        plt.axvline()
                        # plt.savefig('seance2/vp{:.0f}.png'.format(h_tuple[idx]/1000),dpi=120)
                        num_plot+=1

    plt.show()
    plt.close()




# affiche_vp(h_tuple, mac_tuple, km_tuple, ms_tuple,P)




def affiche_vp_faible(h_tuple, mac_tuple, km_tuple, ms_tuple,P):
    num_plot = 1
    fig = plt.figure(1)
    for idx,h in enumerate(h_tuple):
        if num_plot > 8:
            num_plot = 1
            plt.show()
            plt.close()
        for km in km_tuple:
            for ms in ms_tuple:
                for mac in mac_tuple:
                        pt_trim = (h, mac, km, ms)
                        vpp = vp[pt_trim]
                        x = [vpp[k2].real for k2 in range(2, len(vpp))]
                        y = [vpp[k2].imag for k2 in range(2, len(vpp))]
                        # for i,xx in enumerate(x):
                        #         plt.annotate("toto",xy=(xx, y[i]), xytext=(xx-0.2, y[i]-0.2),
                        #                  arrowprops={'facecolor': 'red', 'shrink': 0.05})
                        plt.subplot(4, 2, num_plot)
                        plt.scatter(x, y, color='red',s=4)
                        plt.grid(True)
                        plt.axis([-0.009,0.0,-0.1,0.1])
                        plt.text(-0.008,0.07,"Ma={:.1f}, km = {:.1f} , ms = {:.1f}".format(pt_trim[1],pt_trim[2],pt_trim[3]), fontsize=8)
                        fig.subplots_adjust(hspace=0.3)
                        plt.suptitle("Valeurs propres pour h = {:.0f} km pour le {}".format(h_tuple[idx]/1000,P.name))
                        print("pt_trim = ", pt_trim, " \n valeurs propres ", vpp, "\n")
                        plt.xticks(fontsize=7)
                        plt.yticks(fontsize=7)
                        plt.axhline()
                        plt.axvline()
                        # plt.savefig('seance2/vp{:.0f}_faible.png'.format(h_tuple[idx]/1000),dpi=120)
                        num_plot+=1

    plt.show()
    plt.close()




affiche_vp_faible(h_tuple, mac_tuple, km_tuple, ms_tuple,P)


# Pour la séance 3, le point de trim
# pt_trim =  (11000.0, 0.8, 0.9, 1.0)  a toutes ses vp à partie réelle <0
# valeurs propres  [
# -0.64246345 + 2.71414129j \
# -0.64246345 - 2.71414129j \
# -0.00269042 + 0.05730786j
# -0.00269042 - 0.05730786j
# ]



# pt_trim=(11000.0, 0.8, 0.9, 1.0)
# vpp = vp[pt_trim]
# x = [vpp[i].real for i in range(0,len(vpp))]
# y = [vpp[i].imag for i in range(0,len(vpp))]
# print("real", x)
# print("imag", y)
# plt.scatter(x,y,color='red')
# plt.grid(True)
# plt.title("trim (h, Mach, km, ms) : {}".format(pt_trim))
# plt.show()
