# -*- coding: utf-8 -*-
# /usr/bin/python3.5 /home/blackdoe/Dropbox/TSP/CITI/stage1A/resultatsFinaux/P10.py

# Optimisation de réseau d'antennes: ambiguïtés et SNR critique
# Premier problème préliminaire: détermination du réseau optimal

from numpy import sqrt, pi, cos, linspace, amax, amin, zeros, where, squeeze
import matplotlib.pyplot as plt

Npoints, N = 1000, 100
Ndeltas = 2 * N
w = linspace(0, pi, Npoints, endpoint=True)
d = linspace(-0.15, 0, Ndeltas, endpoint=True)

conversion_indice_delta = lambda a: 0.3 * a / (Ndeltas - 1) - 0.15

# Cas P = 10 ###########################################################################################################

print('\nCas P=10\n')

ap1_u, ap2_u, ap3_u, ap4_u, ap5_u = 1/4, 3/4, 5/4, 7/4, 9/4
ap1 = ap1_u
ap2 = ap2_u
ap3 = ap3_u
ap4 = ap4_u
ap5 = sqrt(165/16 - ap1**2 - ap2**2 - ap3**2 - ap4**2)
g10 = lambda w: 4 * (cos(2*ap1*w) + cos(2*ap2*w) + cos(2*ap3*w) + cos(2*ap4*w) + cos(2*ap5*w)) ** 2
g10_u = lambda w: 4 * (cos(2*ap1_u*w) + cos(2*ap2_u*w) + cos(2*ap3_u*w) + cos(2*ap4_u*w) + cos(2*ap5_u*w)) ** 2

# Détermination du cas optimisé:

A10 = zeros((Ndeltas,Ndeltas,Ndeltas,Ndeltas))  # Matrice de stockage des maxima de lobes secondaires.
# Ligne i : même valeur de delta1
# Colonne j: même valeur de delta2
# Dimensions k,l: mêmes valeur de delta3, delta4 respectivement

i,j,k, l = 0,0,0,0
for delta1 in d:
    for delta2 in d:
        for delta3 in d:
            for delta4 in d:
                ap1 = ap1_u + delta1
                ap2 = ap2_u + delta2
                ap3 = ap3_u + delta3
                ap4 = ap4_u + delta4
                ap5 = sqrt(165/16 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2 - ap4 ** 2)
                u = g10(w)
                u[w < 2 * pi/10] = 0  # Annule les valeurs du lobe principal. Critère approximatif
                A10[i,j,k,l] = amax(u)
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

u = g10_u(w)
u[w < 2 * pi/10] = 0
A10u = amax(u)

# Détermination du delta optimal, calcul des api optimisés en fonction de la position dans A8 du min des amplitudes

delta1_o = conversion_indice_delta(squeeze(where(A10 == amin(A10)))[0])
delta2_o = conversion_indice_delta(squeeze(where(A10 == amin(A10)))[1])
delta3_o = conversion_indice_delta(squeeze(where(A10 == amin(A10)))[2])
delta4_o = conversion_indice_delta(squeeze(where(A10 == amin(A10)))[3])

ap1 = ap1_u + delta1_o
ap2 = ap2_u + delta2_o
ap3 = ap3_u + delta3_o
ap4 = ap4_u + delta4_o
ap5 = sqrt(165/16 - ap1 ** 2 - ap2 ** 2 - ap3 ** 2 - ap4 ** 2)

print('ap1 = ', ap1, '\nap2 = ', ap2, '\nap3 = ', ap3, '\nap4 = ', ap4, '\nap5 = ', ap5)

# Résultat

print('\nRéseau uniforme:\nA_0 = ',(10 ** 2),'\nA_1 = ',amin(A10u),'\nA_0 / A_1 = ', (10 ** 2)/amin(A10u))
print('\nRéseau optimisé:\nA_0 = ',(10 ** 2),'\nA_1 = ',amin(A10),'\nA_0 / A_1 = ', (10 ** 2)/amin(A10))

# Tracé de g10(omega) en log sur y

figP10 = plt.figure()

# Tracé pour le cas uniforme

plt.semilogy(w, g10_u(w), color='black', label='Cas uniforme: $\delta = 0$')

# Tracé pour le cas optimisé

plt.semilogy(w, g10(w), label='Cas optimisé: $\delta_1 = {:.3f}, \delta_2 = {:.3f}, \delta_3 = {:.3f}, \delta_4 = {:.3f}$'.format(delta1_o,delta2_o,delta3_o,delta4_o))

plt.ylim((1e-4, 1e2))
plt.title('Cas $P=10$')
plt.hold(True)
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')

########################################################################################################################

# Tracé final de toutes les courbes

plt.show()
