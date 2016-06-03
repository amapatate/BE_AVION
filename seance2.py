#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import math
import numpy as np, matplotlib.pyplot as plt
import utils as ut
import dynamic as dy

markeur = ('o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', ' ')
plt.grid(True)
plt.axhline()
plt.axvline()
P = dy.Param_737_300()

#
# def plot_trim():
#     plt.title("point de trim en fonction  de l'altitude : " + P.name)
#     # plt.axis([2000., 12000., 0., 30])
#     #####################################################################################
#     # FONCTIONS ANNEXES : MATHPLOTLIB
#     # plt.text(0.52, 70, "La poussée décroît avec l'altitude et le mach")
#     # plt.annotate('3000m', xy=(0.7, 89), xytext=(0.65, 75),
#     #              arrowprops={'facecolor': 'red', 'shrink': 0.05})
#     ####################################################################################
#     # for idx, h in enumerate(h_list):
#     #     X = [0., h, 0., 0., 0., 0.]  # (y h va alpha theta q)
#     #     ord_F = np.array([dynamic.propulsion_model_mach(X, U, P, mach) for mach in abs_mach]) / 1000.  # en kN
#     #     plt.plot(abs_mach, ord_F, marker=markeur[idx], label=str(int(h)) + "m")
#     plt.plot(h_tuple, alphas)
#     plt.legend()
#     plt.xlabel('mach')
#     plt.ylabel('Poussee max en kN')
#     # plt.savefig('images/poussee_of_mach.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
#     plt.show()


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


# for i, km in enumerate(km_tuple):
#     for j, ms in enumerate(ms_tuple):
#         trace(km, ms)
#         plt.show()
#         plt.close()

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
print("angle réglage PHR dphr16 = ", dphr16, " rad" )
print("Gaz ou thrust dth16 = ", dth16)
print("CL16 = ", CL16)



# 4.2.2.4 Conclure quant à la technique de réglage de la vitesse de l’avion.

# >>>> pour augmenter la vitesse, il faut diminuer CL donc alphae donc augmenter algébriquement
# l'angle du PHR pour augmenter la portance ou diminuer la déportance.

# 4.2.2.5 - Utiliser la valeur de CL obtenue ainsi que les graphiques tracés lors de la première séance
# pour déterminer (approximativement) les valeurs de trim αe et δPHRe .

print("La valeur de alphae obtenue via le graphique de la séance 1 CLe=f(alphae) est ", ut.rad_of_deg(6.7))
print("La valeur de dphre obtenue via le graphique dphre = f(alphae) de la séance 1 en utilisant la valeur de alphae précédente est -0.21rad")

# 4.2.2.6 - Comment obtenir la valeur de trim de la manette des gaz δth à partir de la polaire équilibrée ?
# connaissant CLe = 0.57 la polaire équilibrée donne CDe = 0.04  ; on utilise F = D poussee = trainee pour calculé deltath

# calcul du coef de poussée tel que F = coef * deltath
coef = dy.coef_pousseel(X,U,P)

#calcul de la trainée
CDe = 0.04
D = dy.trainee(h,mac, CDe,P)

# deltath
deltath_polaire = D/coef

print("deltath obtenu via la polaire équilibrée = ", deltath_polaire )