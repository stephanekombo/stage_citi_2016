
 - A priori OK pour P = 4,5,6,7,8,9. Voir P456789.py
	Voir resultats.log.txt et fichiers PX.png pour résultats (valeurs des ap "non" arrondies)

 - P = 10: voir P10.py

  * Résultats étonnants (voir P10.png), malgré la précision (200 valeurs des $\delta_i$)
	Voir resultatsP10.log.txt pour valeurs numériques.
  * Pour tenter avec plus de précision, éditer N (ligne 10) et adapter ligne 15:
	conversion_indice_delta est une fonction affine renvoyant,
	pour chaque valeur de [0, Ndeltas - 1], une des Ndelta valeurs decimales entre les bornes
	choisies dans la définition de d (ligne 13)

 - P = 11: voir P11.py

  * Non traité (temps de calcul trop long). Très similaire à P11.py
