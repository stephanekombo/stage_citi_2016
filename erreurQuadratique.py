# -*- coding: utf-8 -*-
# /usr/bin/python3.5 /home/blackdoe/Dropbox/TSP/CITI/stage1A/resultatsFinaux/seuil/erreurQuadratique.py

# Optimisation de réseau d'antennes: ambiguïtés et SNR critique
# Deuxième problème préliminaire: calcul du seuil critique

from numpy import linspace, array, exp, pi, sin, arange, absolute, dot, random, argmax, vectorize, cos, mean, empty, zeros
import matplotlib.pyplot as plt

# Paramétrisation ######################################################################################################

P = 6
I = 100
K1 = 200
K2 = 500
n = 1  # Différentes valeurs de $\sigma_n^2$.
t0 = 0  # t0 = $\theta_0$
T = linspace(-pi/2, pi/2, I, endpoint=True)  # t = $\theta$

s = array([0.01, 0.02, 0.03, 0.06, 0.1, 0.2, 0.3, 0.6, 1, 2, 3, 6, 10, 20, 30, 60, 100])  # Valeurs de $\sigma_s^2$.
# Tatonner selon les valeurs de P. Examiner en particulier [1e-1,1e1]

# Positions des réseaux optimisés. Note: a2Q = [-aQ, ..., -a1, a1, ... aQ], a2Q+1 = [-aQ, ..., -a1,0, a1, ... aQ]
# Arrondies à 1e-4 près. Arbitraire...

ap4 = array([-0.7634, -0.2055, 0.2055, 0.7634])
ap5 = array([-1.0113, -0.4766, 0, 0.4766, 1.0113])
ap6 = array([-1.2868, -0.692, -0.2296, 0.2296, 0.692, 1.2868])
ap7 = array([-1.5466, -0.954, -0.445, 0, 0.445, 0.954, 1.5466])
ap8 = array([-1.8256, -1.1935, -0.6678, -0.2161, 0.2161, 0.6678, 1.1935, 1.8256])
ap9 = array([-2.0991, -1.4495, -0.8952, -0.4374, 0, 0.4374, 0.8952, 1.4495, 2.0991])
ap10 = array([-2.2592, -1.8111, -1.1452, -0.7206, -0.3126, 0.3126, 0.7206, 1.1452, 1.8111, 2.2592])
ap11 = array([-5/2., -2., -3/2., -1., -1/2., 0.,  1/2., 1., 3/2., 2., 5/2.])  # Réseau uniforme en absence de résultats

ap = {4: ap4, 5: ap5, 6: ap6, 7: ap7, 8: ap8, 9: ap9, 10: ap10, 11: ap11}  # Pour utiliser ap[P] dans la suite

# Position des réseaux uniformes

ap4u = array([-3/4., -1/4., 1/4., 3/4.])
ap5u = array([-1., -1/2., 0., 1/2., 1.])
ap6u = array([-5/4., -3/4., -1/4., 1/4., 3/4., 5/4.])
ap7u = array([-3/2., -1., -1/2., 0., 1/2., 1., 3/2.])
ap8u = array([-7/4., -5/4., -3/4., -1/4., 1/4., 3/4., 5/4., 7/4.])
ap9u = array([-2., -3/2., -1., -1/2., 0.,  1/2., 1., 3/2., 2.])
ap10u = array([-9/4., -7/4., -5/4., -3/4., -1/4., 1/4., 3/4., 5/4., 7/4., 9/4.])
ap11u = array([-5/2., -2., -3/2., -1., -1/2., 0.,  1/2., 1., 3/2., 2., 5/2.])

apu = {4: ap4u, 5: ap5u, 6: ap6u, 7: ap7u, 8: ap8u, 9: ap9u, 10: ap10u, 11: ap11u}

# Autres fonctions utiles

ah = lambda t,P: exp(2j * pi * sin(-t) * ap[P])
ahu = lambda t,P: exp(2j * pi * sin(-t) * apu[P])
L = lambda X,i,t,K: mean(array([absolute(dot(ah(t,P), X[i,k]))**2 for k in arange(K)]))  # Renvoie $K.L_x(t=\theta)$
Lu = lambda X,i,t,K: mean(array([absolute(dot(ahu(t,P), X[i,k]))**2 for k in arange(K)]))  # Idem, cas uniforme

# Calcul de l'erreur quadratique #######################################################################################

def e(s,ah,L,P,K):

    T_LM = empty(I)
    S = empty((I,K),dtype=complex)
    N = X = empty((I,K,P),dtype=complex)  # I simulations de K vecteurs aléatoires $(X_k, k=1..K)$ de taille P

    for i in arange(I):
        for k in arange(K):
            S[i,k] = random.normal(0, s/2) + random.normal(0, s/2) * 1j
            N[i,k] = random.normal(0, n/2, P) + random.normal(0, n/2, P) * 1j
            X[i,k] = S[i,k] * ah(t0,P) + N[i,k]
        T_LM[i] = argmax(array([L(X,i,t,K) for t in T]))

    T_LM = T_LM * pi / (I - 1) - pi / 2  # Convertit les indices argmax en les valeurs de $\theta_i$ correspondantes

    return mean((T_LM - t0)**2)

E = vectorize(e) # Règle le problème de types/dimensions

# Définition des moments pour chaque valeur de P et de la borne de Cramer-Rao

Sx = {4: 5/4, 5: 5/2, 6: 35/8, 7: 7, 8: 21/2, 9: 15, 10: 165/8, 11: 55/2}
f = lambda x,K,P: (1 + x * P) * ( 2*K*(x**2)*P*4*(pi**2)*(cos(t0)**2)*Sx[P]) ** (-1)  # x = s/n

# Tracé des courbes en log/log #########################################################################################

# Cas uniforme et optimisé de front, avec tracé de la borne de Cramer-Rao

fig, (u, o) = plt.subplots(1, 2, sharex=True, sharey=True)

u.loglog(s/n, E(s, ahu, Lu, P, K1), label='Cas uniforme, K1', linestyle='--', marker='.')  # Try marker='x'
u.loglog(s/n, E(s, ahu, Lu, P, K2), label='Cas uniforme, K2', linestyle='--', marker='.')
o.loglog(s/n,E(s,ah,L,P,K1), label='Cas optimisé, K1', linestyle='--', marker='.')
o.loglog(s / n, E(s, ah, L, P, K2), label='Cas optimisé, K2', linestyle='--', marker='.')

u.loglog(s/n,f(s/n,K1,P), color='black', label='Borne de Cramer Rao, K1', linestyle='solid')
u.loglog(s/n,f(s/n,K2,P), color='brown', label='Borne de Cramer Rao, K2', linestyle='solid')
o.loglog(s/n,f(s/n,K1,P), color='black', label='Borne de Cramer Rao, K1', linestyle='solid')
o.loglog(s/n,f(s/n,K2,P), color='brown', label='Borne de Cramer Rao, K2', linestyle='solid')

#plt.xlim((1e-10,1e2))
#plt.ylim((1.5*1e-1, 4.5*1e3))
u.set_title('Cas uniforme: $P = {:.0f}, I = {:.0f}, K_1 = {:.0f}, K_2 = {:.0f}$'.format(P,I,K1,K2))
o.set_title('Cas optimisé: $P = {:.0f}, I = {:.0f}, K_1 = {:.0f}, K_2 = {:.0f}$'.format(P,I,K1,K2))
u.hold(True)
o.hold(True)
u.grid(True)
o.grid(True)
u.legend(loc='lower left', fontsize='small')
o.legend(loc='lower left', fontsize='small')

# Tracé final de toutes les courbes

plt.show()
