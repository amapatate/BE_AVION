#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np, matplotlib.pyplot as plt
import utils, dynamic

markeur = ('o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd')
plt.grid(True)
P = dynamic.Param_737_300()


# 3. Séance 1
# 3.1 Travail hors séance encadrée : 2 diapos sur le 737_300
# année de mise en service, nombre d'exemplaires construits, compagnies exploitantes,capacité, performances, etc

# 3.2 Séance encadrée
# 3.1.2 Tracer la poussée max F(Mach) avec Mach entre 0.5 et 0.8 pour h = 3000m et 11000m
# comment évolue la poussée max avec l'altitude et le nombre de mach

def poussee():
    abs_mach = np.linspace(0.5, 0.8, 20)  # mach en abscissse : vecteur avec 20 points
    U = [0., 1.]  # [ delta_phr delta_thrust ]  poussé max donc thrust=1
    h_list = [3000., 10000.]
    plt.title("Poussée max en fonction du mach et de l'altitude : "+P.name)
    plt.axis([0.5, 0.8, 50, 100])
    #####################################################################################
    # FONCTIONS ANNEXES : MATHPLOTLIB
    plt.text(0.52, 70, "La poussée décroît avec l'altitude et le mach")
    plt.annotate('3000m', xy=(0.7, 89), xytext=(0.65, 75),
                 arrowprops={'facecolor': 'red', 'shrink': 0.05})
    ####################################################################################
    for idx, h in enumerate(h_list):
        X = [0., h, 0., 0., 0., 0.]  # (y h va alpha theta q)
        ord_F = np.array([dynamic.propulsion_model_mach(X, U, P, mach) for mach in abs_mach]) / 1000.  # en kN
        plt.plot(abs_mach, ord_F, marker=markeur[idx], label=str(int(h)) + "m")

    plt.legend()
    plt.xlabel('mach')
    plt.ylabel('Poussee max en kN')
    plt.savefig('images/poussee_of_mach.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()


#poussee()


# 3.2.2 Tracer le coefficient de portance CL en fonction de l’incidence α (comprise entre −10 π/180 rad
# à +20 π/180 rad ) lorsque δPHR vaut −30pi/180 rad puis +20pi/180 rad. Quel est l’effet de δPHR
# sur le coefficient de portance ? Le modèle proposé simule t-il le décrochage de l’avion ?
def portance():

    dphr_deg=[-30.,20.]


    dphr_rad= np.array(dphr_deg) * np.pi/180.

    a_deg=[-10,20]
    a_rad= np.array(a_deg) * np.pi/180.

    abs_alpha = np.linspace(a_rad[0], a_rad[1], 20)  # mach en abscissse : vecteur avec 20 points
    abs_alpha_deg = np.linspace(a_deg[0], a_deg[1], 20)
    plt.title("Coefficient de portance : "+P.name)
    plt.axis([-20,30,-2,3])
    #####################################################################################
    # FONCTIONS ANNEXES : MATHPLOTLIB
    plt.text(-19.5, 2.6, "Décrochage non pris en compte par le modèle")
    plt.text(-19.5, 2.4, "La Portance croît avec le phr.")
    plt.annotate("Modèle linéaire : Pas de chute de portance : $\delta$phr = " + str(int(dphr_deg[1]))+"°",
    xy=(20, 2.8), xytext=(-19.5, 2.8),arrowprops={'facecolor':'red', 'shrink':0.05} )
    ####################################################################################
    for idx, dphr in enumerate(dphr_rad):
        ord_CL = np.array([dynamic.get_aero_ceofs(100., alpha, 0., dphr, P)[0] for alpha in abs_alpha])
        plt.plot(abs_alpha_deg, ord_CL, marker=markeur[idx], label="$\delta$phr = "+ str(int(dphr_deg[idx]))+ "°")

    plt.legend(loc=4)
    plt.axhline()
    plt.axvline()
    plt.xlabel(r'Incidence $\alpha$')
    plt.ylabel('Coef. portance CL')
    plt.savefig('images/CL_of_alpha.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()


#portance()



# 3.2.3
# Tracer le coefficient Cm du moment de tangage en fonction de α pour quatre valeurs de
# la marge statique ms : −0.1, 0, 0.2 et 1 lorsque δPHR = 0. Quel est la réaction de l’avion
# en cas d’augmentation intempestive de l’incidence α ?


def tangage():

    a_deg=[-10,20]
    a_rad= np.array(a_deg) * np.pi/180.
    ms_list = [-0.1,0.,0.2,1]
    abs_alpha = np.linspace(a_rad[0], a_rad[1], 20)  # mach en abscissse : vecteur avec 20 points
    abs_alpha_deg = np.linspace(a_deg[0], a_deg[1], 20)
    plt.title("Coef. de tangage pour diverses marges statiques :"+P.name)
    # plt.axis([-20,30,-2,3])
    #####################################################################################
    # FONCTIONS ANNEXES : MATHPLOTLIB
    plt.text(-9, -1.5, "ms = 1, Cm décroît donc avion stable statiquement")
    plt.text(-9, -1.9, "ms = -0.1, Cm croît donc avion instable statiquement")
    # plt.text(-19.5, 2.4, "La Portance croît avec le phr.")
    # plt.annotate("Modèle linéaire : Pas de chute de portance : $\delta$phr = " + str(int(dphr_deg[1]))+"°",
    # xy=(20, 2.8), xytext=(-19.5, 2.8),arrowprops={'facecolor':'red', 'shrink':0.05} )
    ####################################################################################
    for idx, ms in enumerate(ms_list):
        ord_CL = np.array([dynamic.get_aero_ceofs_ms(100., alpha, 0., 0., P, ms)[2] for alpha in abs_alpha])
        plt.plot(abs_alpha_deg, ord_CL, marker=markeur[idx], label="ms = "+ str(ms_list[idx]))

    plt.legend(loc=4)
    plt.axhline()
    plt.axvline()
    plt.xlabel(r'Incidence $\alpha$')
    plt.ylabel('Coef. tangage Cm')
    plt.savefig('images/Cm_of_alpha.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()

#tangage()


#3.2.4
# Calculer et tracer en fonction de l’incidence α la valeur δPHRe de δPHR pour laquelle le
# moment de tangage est nul (i.e. Cm = 0) et le vol stabilisé (i.e. q = 0). Comment varie
# δPHRe avec la marge statique ms et le volume d’empennage V t ? Comment varie δPHRe
# en fonction de l’incidence d’équilibre que l’on notera αe ?

def dphre(alpha,P):
    return (P.ms * P.CLa * (alpha - P.a0) -P.Cm0)/ P.Cmd

def dphre_ms(alpha,P,mms):
    return (mms * P.CLa * (alpha - P.a0) -P.Cm0)/ P.Cmd

def dphre_of_alpha():
    a_deg=[-10,20]
    a_rad= np.array(a_deg) * np.pi/180.
    ms_list = [-0.1,0.,0.2,1]
    abs_alpha = np.linspace(a_rad[0], a_rad[1], 20)  # mach en abscissse : vecteur avec 20 points
    abs_alpha_deg = np.linspace(a_deg[0], a_deg[1], 20)
    plt.title("$\delta$PHRe à l'équilibre : "+P.name)
    # plt.axis([-20,30,-2,3])
    #####################################################################################
    # FONCTIONS ANNEXES : MATHPLOTLIB
    plt.text(-9, -0.17, r" $\delta$PHRe décroît plus rapidement avec $\alpha$ lorsque ms > 0 croît ")
    plt.text(-9, -0.19, "$\delta$PHRe décroît avec Vt, vol. Empennage")
    plt.text(-9, 0.21, "La Portance croît avec le phr.")
    # plt.annotate("Modèle linéaire : Pas de chute de portance : $\delta$phr = " + str(int(dphr_deg[1]))+"°",
    # xy=(20, 2.8), xytext=(-19.5, 2.8),arrowprops={'facecolor':'red', 'shrink':0.05} )
    ####################################################################################
    for idx, ms in enumerate(ms_list):
        ord_dphre = np.array([dphre_ms(alpha,P, ms) for alpha in abs_alpha])
        plt.plot(abs_alpha_deg, ord_dphre, marker=markeur[idx], label="ms = "+ str(ms_list[idx]))

    plt.legend(loc=3)
    plt.axhline()
    plt.axvline()
    plt.xlabel(r'Incidence $\alpha$')
    plt.ylabel('$\delta$PHRe ')
    plt.savefig('images/dphre_of_alpha.png', dpi=120)  # sauvegarde du graphe au format png dans le dossier images
    plt.show()

dphre_of_alpha()