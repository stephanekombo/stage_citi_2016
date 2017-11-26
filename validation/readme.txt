Voir erreurQuadratique.py

 - Seuls résultats pour des valeurs de I et K trop faibles.
	Il faut relancer les simulations, maintenant que le code a été optimisé.

 - Les images 4.a.png à 4.i.png montrent, avec des zooms successifs, la zone où se trouve le seuil.

  * Il faut relancer la simulation avec le plus de valeurs de s possibles
	comprises entre 0..1 et 10. Attention, très gros temps de calcul à prévoir.
	10 valeurs de s (ligne 20) par simulation me semble convenable.
	Si possible faire tourner pendant le week-end sur un serveur le plus puissant possible.

  * Pas assez de points: borne de Cramer-Rao dépassée...

 - 6..a.png 

  * Résultats moins choquants: Borne de Cramer-Rao respectée.

  * Le seuil semble cette fois compris entre 0.1 et 1.

 - P = 11: voir P11.py

  * Non traité (temps de calcul trop long). Très similaire à P11.py

 - Pour les nouvelles simulations:

  * Paramétrer en modifiant les lignes 12 à 20 du fichier.

  * Selon la précision requise, je recommande FORTEMENT de réduire le plus possible
	le nombre de décimales des positions optimisées (lignes 26 à 33).
	Rappel: le cas P = 10 est suspect

  * Idem, faire tourner un week-end sur un serveur.
