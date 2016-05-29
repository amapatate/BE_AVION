#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np, matplotlib.pyplot as plt
import utils, dynamic

# si le code interprété est dans le fichier courant tuto.py
# alors la vraible spéciale __name__ prend la valeur '__main__'
# si on utilise tuto.py par importation dans un autre fichier __name__ prend la valeur 'tuto.py'
# donc différente de main et le code dans le if ci-dessous n'est pas executé.
if __name__ == "__main__":
    v_Kt = np.array([100, 150, 200])
    v_ms = utils.mps_of_kt(v_Kt)
    print('{} Kt -> {} m/s'.format(v_Kt, v_ms))

    h = 5000
    p, rho, T = utils.isa(h)
    print('atmosphère normalisée à {} m: {:.1f} Pa {:.3} kg/m3 {:.2f} K'.format(h, p, rho, T))

    hs = np.linspace(0, 10000, 21) # tableau numpy de 21 valeurs de 0 à 10000 m [0, 500, 1000, ...
    # isa = np.array([utils.isa(h) for h in hs])

    # isa(h) renvoie un tuple (p,rho,T) isa(h)[0] sélectionne donc la pression p en pascal
    # isa est un tableau numpy ; on divise par 100 pour avoir des hpa
    isa = np.array([utils.isa(h)[0] for h in hs]) / 100
    print(isa)
    plt.plot(hs, isa) # la fonction plot de matplotlib prend le vecteur numpy hs comme coordonnées en abscisse
    # et le vecteur "isa" pour les valeurs correspodantes en ordonnée
    plt.xlabel('altitude en m');
    plt.ylabel('pression en Hectopascal');
    plt.savefig('isa_pressure.png', dpi=120) # sauvegarde du graphe au format png dans le dossier images

    p = dynamic.Param_A320()
    print('longueur du fuselage de l\'{}: {}m'.format(p.name, p.l_fus))

    plt.show()
