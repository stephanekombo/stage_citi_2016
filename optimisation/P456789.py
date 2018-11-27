# -*- coding: utf-8 -*-
# test
# /usr/bin/python3.5 /home/blackdoe/Dropbox/TSP/CITI/stage1A/essaisPreliminaires/resultats/P456789.py

# Optimisation de réseau d'antennes: ambiguïtés et SNR critique
# Premier problème préliminaire: détermination du réseau optimal

import matplotlib.pyplot as plt
from numpy import sqrt, pi, cos, linspace, amax, amin, zeros, where, squeeze

Npoints, N = 1000, 100
Ndeltas = 2 * N
w = linspace(0, pi, Npoints, endpoint=True)
d = linspace(-0.15, 0.15, Ndeltas, endpoint=True)


# Cas P = 4 ############################################################################################################

print('\nCas P=4\n')

ap1_u, ap2_u = 1/4, 3/4
ap1 = ap1_u
ap2 = sqrt(5/8 - ap1 ** 2)
g4 = lambda w: 4 * (cos(2 * ap1 * w) + cos(2 * ap2 * w)) ** 2
g4_u = lambda w: 4 * (cos(2 * ap1_u * w) + cos(2 * ap2_u * w)) ** 2

# Détermination du cas optimisé:

A4 = zeros((Ndeltas))  # Matrice de stockage des maxima de lobes secondaires. Ndeltas valeurs de delta

i = 0
for delta in d:
    ap1 = ap1_u + delta
    ap2 = sqrt(5/8 - ap1 ** 2)
    u = g4(w)
    u[w < 2 * pi/4] = 0  # Annule les valeurs du lobe principal. Critère approximatif à vérifier sur les courbes
    A4[i] = amax(u)
    i += 1

# Détermination du maximum dans le cas uniforme

u = g4_u(w)
u[w < 2 * pi/4] = 0
A4u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A4 du min des amplitudes

conversion_indice_delta = lambda a: 0.3 * a / (Ndeltas - 1) - 0.15

delta_o = conversion_indice_delta(squeeze(where(A4 == amin(A4))))

ap1 = ap1_u + delta_o
ap2 = sqrt(5/8 - ap1 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(4 ** 2),'\nA_1 = ',amin(A4u),'\nA_0 / A_1 = ', (4 ** 2)/amin(A4u))
print('\nRéseau optimisé:\nA_0 = ',(4 ** 2),'\nA_1 = ',amin(A4),'\nA_0 / A_1 = ', (4 ** 2)/amin(A4))

# Tracé de g4(omega) en log sur y

figP4 = plt.figure()

# Tracé pour le cas uniforme

plt.semilogy(w, g4_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g4(w), label='Cas optimisé: $\delta = {:.3f}$'.format(delta_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=4$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Cas P = 5 ############################################################################################################

print('\nCas P=5\n')

ap1_u, ap2_u = 1/2, 1
ap1 = ap1_u
ap2 = sqrt(5/4 - ap1 ** 2)
g5 = lambda w: (1 + 2 * cos(2 * ap1 * w) + 2 * cos(2 * ap2 * w)) ** 2
g5_u = lambda w: (1 + 2 * cos(2 * ap1_u * w) + 2 * cos(2 * ap2_u * w)) ** 2

# Détermination du cas optimisé:

A5 = zeros((Ndeltas))  # Matrice de stockage des maxima de lobes secondaires. Ndeltas valeurs de delta

i = 0
for delta in d:
    ap1 = ap1_u + delta
    ap2 = sqrt(5/4 - ap1 ** 2)
    u = g5(w)
    u[w < 2 * pi/5] = 0  # Annule les valeurs du lobe principal. Critère approximatif à vérifier sur les courbes
    A5[i] = amax(u)
    i += 1

# Détermination du maximum dans le cas uniforme

u = g5_u(w)
u[w < 2 * pi/5] = 0
A5u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A5 du min des amplitudes

delta_o = conversion_indice_delta(squeeze(where(A5 == amin(A5))))

ap1 = ap1_u + delta_o
ap2 = sqrt(5/4 - ap1 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(5 ** 2),'\nA_1 = ',amin(A5u),'\nA_0 / A_1 = ', (5 ** 2)/amin(A5u))
print('\nRéseau optimisé:\nA_0 = ',(5 ** 2),'\nA_1 = ',amin(A5),'\nA_0 / A_1 = ', (5 ** 2)/amin(A5))

# Tracé de g5(omega) en log sur y

figP5 = plt.figure()

# Tracé pour le cas uniforme: ap1 = 1/2; ap2 = 1

plt.semilogy(w, g5_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé: ap1 = ap1_u + delta_o, ap2 = sqrt(5 / 4 - ap1 ** 2)

plt.semilogy(w, g5(w), label='Cas optimisé: $\delta = {:.3f}$'.format(delta_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=5$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Cas P = 6 ############################################################################################################

print('\nCas P=6\n')

ap1_u, ap2_u, ap3_u = 1/4, 3/4, 5/4
ap1 = ap1_u
ap2 = ap2_u
ap3 = sqrt(35/16 - ap1 ** 2 - ap2 ** 2)
g6 = lambda w: 4 * (cos(2 * ap1 * w) + cos(2 * ap2 * w) + cos(2 * ap3 * w)) ** 2
g6_u = lambda w: 4 * (cos(2 * ap1_u * w) + cos(2 * ap2_u * w) + cos(2 * ap3_u * w)) ** 2

# Détermination du cas optimisé:

A6 = zeros((Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires. Ndeltas valeurs de chaque delta i
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2

i,j = 0,0
for delta1 in d:
    for delta2 in d:
        ap1 = ap1_u + delta1
        ap2 = ap2_u + delta2
        ap3 = sqrt(35/16 - ap1 ** 2 - ap2 ** 2)
        u = g6(w)
        u[w < 2 * pi/6] = 0  # Annule les valeurs du lobe principal. Critère approximatif
        A6[i,j] = amax(u)
        j += 1
    i += 1
    j = 0

# Détermination du maximum dans le cas uniforme

u = g6_u(w)
u[w < 2 * pi/6] = 0
A6u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A6 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A6 == amin(A6)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A6 == amin(A6)))[1])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = sqrt(35/16 - ap1 ** 2 - ap2 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(6 ** 2),'\nA_1 = ',amin(A6u),'\nA_0 / A_1 = ', (6 ** 2)/amin(A6u))
print('\nRéseau optimisé:\nA_0 = ',(6 ** 2),'\nA_1 = ',amin(A6),'\nA_0 / A_1 = ', (6 ** 2)/amin(A6))

# Tracé de g6(omega) en log sur y

figP6 = plt.figure()

# Tracé pour le cas uniforme

plt.semilogy(w, g6_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g6(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}$'.format(delta1_o,delta2_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=6$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Cas P = 7 ############################################################################################################

print('\nCas P=7\n')

ap1_u, ap2_u, ap3_u = 1/2, 1, 3/2
ap1 = ap1_u
ap2 = ap2_u
ap3 = sqrt(7/2 - ap1 ** 2 - ap2 ** 2)
g7 = lambda w: (1 + 2 * cos(2 * ap1 * w) + 2 * cos(2 * ap2 * w) + 2 * cos(2 * ap3 * w)) ** 2
g7_u = lambda w: (1 + 2 * cos(2 * ap1_u * w) + 2 * cos(2 * ap2_u * w) + 2 * cos(2 * ap3_u * w)) ** 2

# Détermination du cas optimisé:

A7 = zeros((Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires. Ndeltas valeurs de chaque delta i
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2

i,j = 0,0
for delta1 in d:
    for delta2 in d:
        ap1 = ap1_u + delta1
        ap2 = ap2_u + delta2
        ap3 = sqrt(7/2 - ap1 ** 2 - ap2 ** 2)
        u = g7(w)
        u[w < 2 * pi/7] = 0  # Annule les valeurs du lobe principal. Critère approximatif
        A7[i,j] = amax(u)
        j += 1
    i += 1
    j = 0

# Détermination du maximum dans le cas uniforme

u = g7_u(w)
u[w < 2 * pi/7] = 0
A7u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A7 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A7 == amin(A7)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A7 == amin(A7)))[1])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = sqrt(7/2 - ap1 ** 2 - ap2 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(7 ** 2),'\nA_1 = ',amin(A7u),'\nA_0 / A_1 = ', (7 ** 2)/amin(A7u))
print('\nRéseau optimisé:\nA_0 = ',(7 ** 2),'\nA_1 = ',amin(A7),'\nA_0 / A_1 = ', (7 ** 2)/amin(A7))

# Tracé de g7(omega) en log sur y

figP7 = plt.figure()

# Tracé pour cas uniforme

plt.semilogy(w, g7_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g7(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}$'.format(delta1_o,delta2_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=7$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Cas P = 8 ############################################################################################################

print('\nCas P=8\n')

ap1_u, ap2_u, ap3_u, ap4_u = 1/4, 3/4, 5/4, 7/4
ap1 = ap1_u
ap2 = ap2_u
ap3 = ap3_u
ap4 = sqrt(63/12 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)
g8 = lambda w: 4 * (cos(2 * ap1 * w) + cos(2 * ap2 * w) + cos(2 * ap3 * w) + cos(2 * ap4 * w)) ** 2
g8_u = lambda w: 4 * (cos(2 * ap1_u * w) + cos(2 * ap2_u * w) + cos(2 * ap3_u * w) + cos(2 * ap4_u * w)) ** 2

# Détermination du cas optimisé:

A8 = zeros((Ndeltas,Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires.
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2
# Dimension k: même valeur de delta3

i,j,k = 0,0,0
for delta1 in d:
    for delta2 in d:
        for delta3 in d:
            ap1 = ap1_u + delta1
            ap2 = ap2_u + delta2
            ap3 = ap3_u + delta3
            ap4 = sqrt(63/12 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)
            u = g8(w)
            u[w < 2 * pi/8] = 0  # Annule les valeurs du lobe principal. Critère approximatif
            A8[i,j,k] = amax(u)
            k += 1
        j += 1
        k = 0
    i += 1
    j = 0
    k = 0

# Détermination du maximum dans le cas uniforme

u = g8_u(w)
u[w < 2 * pi/8] = 0
A8u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A8 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A8 == amin(A8)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A8 == amin(A8)))[1])
delta3_o = conversion_indice_delta(squeeze(where(A8 == amin(A8)))[2])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = ap3_u + delta3_o
ap4 = sqrt(63/12 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3, '\nap4 = ', ap4)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(8 ** 2),'\nA_1 = ',amin(A8u),'\nA_0 / A_1 = ', (8 ** 2)/amin(A8u))
print('\nRéseau optimisé:\nA_0 = ',(8 ** 2),'\nA_1 = ',amin(A8),'\nA_0 / A_1 = ', (8 ** 2)/amin(A8))

# Tracé de g8(omega) en log sur y

figP8 = plt.figure()

# Tracé pour le cas uniforme

plt.semilogy(w, g8_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g8(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}, \delta_3 = {:.3f}$'.format(delta1_o,delta2_o,delta3_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=8$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Cas P = 9 ############################################################################################################

print('\nCas P=9\n')

ap1_u, ap2_u, ap3_u, ap4_u = 1/2, 1, 3/2, 2
ap1 = ap1_u
ap2 = ap2_u
ap3 = ap3_u
ap4 = sqrt(15/2 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)
g9 = lambda w: (1 + 2 * cos(2 * ap1 * w) + 2 * cos(2 * ap2 * w) + 2 * cos(2 * ap3 * w) + 2 * cos(2 * ap4 * w)) ** 2
g9_u = lambda w: (1 + 2 * cos(2*ap1_u*w) + 2 * cos(2*ap2_u*w) + 2 * cos(2*ap3_u*w) + 2 * cos(2*ap4_u*w)) ** 2

# Détermination du cas optimisé:

A9 = zeros((Ndeltas,Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2
# Dimension k: même valeur de delta3

i,j,k = 0,0,0
for delta1 in d:
    for delta2 in d:
        for delta3 in d:
            ap1 = ap1_u + delta1
            ap2 = ap2_u + delta2
            ap3 = ap3_u + delta3
            ap4 = sqrt(15/2 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)
            u = g9(w)
            u[w < 2 * pi/9] = 0  # Annule les valeurs du lobe principal. Critère approximatif
            A9[i,j,k] = amax(u)
            k += 1
        j += 1
        k = 0
    i += 1
    j = 0
    k = 0

# Détermination du maximum dans le cas uniforme

u = g9_u(w)
u[w < 2 * pi/9] = 0
A9u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A9 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A9 == amin(A9)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A9 == amin(A9)))[1])
delta3_o = conversion_indice_delta(squeeze(where(A9 == amin(A9)))[2])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = ap3_u + delta3_o
ap4 = sqrt(15/2 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3, '\nap4 = ', ap4)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(9 ** 2),'\nA_1 = ',amin(A9u),'\nA_0 / A_1 = ', (9 ** 2)/amin(A9u))
print('\nRéseau optimisé:\nA_0 = ',(9 ** 2),'\nA_1 = ',amin(A9),'\nA_0 / A_1 = ', (9 ** 2)/amin(A9))

# Tracé de g9(omega) en log sur y

figP9 = plt.figure()

# Tracé pour cas uniforme

plt.semilogy(w, g9_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g9(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}, \delta_3 = {:.3f}$'.format(delta1_o,delta2_o,delta3_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=9$')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


# Tracé final de toutes les courbes

plt.show()