# -*- coding: utf-8 -*-
# /usr/bin/python3.5 /home/blackdoe/Dropbox/TSP/CITI/stage1A/resultatsFinaux/P11.py

# Optimisation de réseau d'antennes: ambiguïtés et SNR critique
# Premier problème préliminaire: détermination du réseau optimal

from numpy import sqrt, pi, cos, linspace, amax, amin, zeros, where, squeeze
import matplotlib.pyplot as plt

Npoints, N = 1000, 100
Ndeltas = 2 * N
w = linspace(0, pi, Npoints, endpoint=True)
d = linspace(-0.15, 0, Ndeltas, endpoint=True)

conversion_indice_delta = lambda a: 0.3 * a / (Ndeltas - 1) - 0.15

# Cas P = 11 ###########################################################################################################

print('\nCas P=11\n')

ap1_u, ap2_u, ap3_u, ap4_u, ap5_u = 1/2, 1, 3/2, 2, 5/2
ap1 = ap1_u
ap2 = ap2_u
ap3 = ap3_u
ap4 = ap4_u
ap5 = sqrt(55/4 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2 - ap4 ** 2)
g11 = lambda w: (1 + 2 * cos(2*ap1*w) + 2 * cos(2*ap2*w) + 2 * cos(2*ap3*w) + 2 * cos(2*ap4*w) + 2 * cos(2*ap5* w)) ** 2
g11_u = lambda w: (1+ 2*cos(2*ap1_u*w) + 2*cos(2*ap2_u*w) + 2*cos(2*ap3_u*w) + 2*cos(2*ap4_u*w) + 2*cos(2*ap5_u*w)) ** 2

# Détermination du cas optimisé:

A11 = zeros((Ndeltas,Ndeltas,Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires.
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2
# Dimensions k,l: mêmes valeur de delta3, delta4 respectivement

i,j,k,l = 0,0,0,0
for delta1 in d:
    for delta2 in d:
        for delta3 in d:
            for delta4 in d:
                ap1 = ap1_u + delta1
                ap2 = ap2_u + delta2
                ap3 = ap3_u + delta3
                ap4 = ap4_u + delta4
                ap5 = sqrt(55/4 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2 - ap4 ** 2)
                u = g11(w)
                u[w < 2 * pi/11] = 0  # Annule les valeurs du lobe principal. Critère approximatif
                A11[i,j,k,l] = amax(u)
                l += 1
            k += 1
            l = 0
        j += 1
        k = 0
        l = 0
    i += 1
    j = 0
    k = 0
    l = 0

# Détermination du maximum dans le cas uniforme

u = g11_u(w)
u[w < 2 * pi/11] = 0
A11u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A11 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A11 == amin(A11)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A11 == amin(A11)))[1])
delta3_o = conversion_indice_delta(squeeze(where(A11 == amin(A11)))[2])
delta4_o = conversion_indice_delta(squeeze(where(A11 == amin(A11)))[3])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = ap3_u + delta3_o
ap4 = ap4_u + delta4_o
ap5 = sqrt(55/4 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2 - ap4 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3, '\nap4 = ', ap4, '\nap5 = ', ap5)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(11 ** 2),'\nA_1 = ',amin(A11u),'\nA_0 / A_1 = ', (11 ** 2)/amin(A11u))
print('\nRéseau optimisé:\nA_0 = ',(11 ** 2),'\nA_1 = ',amin(A11),'\nA_0 / A_1 = ', (11 ** 2)/amin(A11))

# Tracé de g11(omega) en log sur y

figP11 = plt.figure()

# Tracé pour le cas uniforme

plt.semilogy(w, g11_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g11(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}, \delta_3 = {:.3f}, \delta_4 = {:.3f}$'.format(delta1_o,delta2_o,delta3_o,delta4_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=11$')
plt.hold(True)
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')

########################################################################################################################

# Tracé final de toutes les courbes

plt.show()
