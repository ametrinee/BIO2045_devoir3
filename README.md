# Stratégie de vaccination spatiale : simulation d'une intervention épidémique sous contraintes

## Organisation du projet

Ce projet a pour objectif de simuler la propagation spatiale d'une maladie infectieuse mortelle et d'évaluer l'efficacité d'une campagne de vaccination ciblée. Il utilise un modèle épidémiologique spatialisé basé sur des agents (individus) se déplaçant sur une grille de coordonnées. L'intervention évaluée repose sur une stratégie qui utilise des tests antigéniques rapides (RAT) pour détecter les foyers d'infection asymptomatiques et construire des barrières immunitaires spatiales.

## Objectifs du projet

L'objectif principal est de maximiser la réduction de la mortalité au sein d'une population initialement naïve de 3750 individus, répartie sur une grille spatiale délimitée (-50 à 50 sur les axes X et Y). La difficulté majeure de cette simulation réside dans la nature insidieuse de la maladie : elle est strictement asymptomatique, dure 21 jours et s'avère 100 % fatale, avec une probabilité de transmission de 0.4 lors d'un contact. À cette dynamique épidémique s'ajoutent des contraintes logistiques et financières rigoureuses. Disposant d'un budget strictement limité à 21 000 $, l'intervention doit optimiser l'utilisation de tests antigéniques rapides (RAT) à 4 $ l'unité (offrant 95 % de précision) et de vaccins à 17 $, dont l'immunité n'est acquise que deux jours après l'inoculation. De plus, aucune mesure sanitaire ne peut être déployée avant le décès du premier cas index. Dans ce contexte de ressources restreintes et d'information imparfaite, le projet s'attache à évaluer d'une part l'efficacité de la stratégie de vaccination mise en place (via la diminution de la mortalité par rapport à un scénario sans intervention), et d'autre part le coût total de cette campagne ciblée.

## Modèle utilisé

Le modèle se base sur des agents stochastiques se déplaçant aléatoirement sur la grille. La transmission ne se produit que lorsque deux agents se trouvent dans la même cellule. La stratégie simulée consiste à administrer des tests RAT à un sous-ensemble de la population et, en cas de résultat positif, à vacciner immédiatement les individus se trouvant dans la même cellule (voisinage immédiat) pour bloquer la chaîne de transmission locale.

## Implémentation

L'implémentation du modèle repose sur le suivi individuel de chaque agent, intégrant son état de santé, son statut vaccinal et son historique de dépistage. La simulation permet de :

- Simuler la dynamique spatio-temporelle de l'épidémie.
- Gérer la décroissance du budget à chaque test ou vaccin administré.
- Répliquer les simulations de manière stochastique pour capturer la variabilité des résultats.
- Générer des graphiques (séries temporelles) illustrant la dynamique de la population (Sains, Infectieux, Décédés).

## Résultat

Cette section présente les données obtenues après 100 réplications stochastiques de la simulation (Population initiale N=3750, Budget = 21 000 $).

| Indicateur | Valeur moyenne (± Écart-type) |
| :--- | :--- |
| **Survivants (S)** | **1 692** (± **805**) |
| **Décès (D)** | **2 058** (± **805**) |
| **Budget restant** | **3 101 $** (± **5 317 $**) |
| **Statut final** | **Extinction complète (I = 0)** |

En moyenne, 1 692 individus ont survécu à l'épidémie, ce qui représente un taux de survie de 45,1 % (écart-type : 805). Sur le plan financier, la campagne de vaccination en anneau a consommé en moyenne 17 899&nbsp;$, laissant un solde budgétaire final de 3 101&nbsp;$ (écart-type : 5 317$). Dans l'ensemble des 100 simulations, la chaîne de transmission s'est interrompue avant l'épuisement total de la population saine.

## Conclusion

Ce projet démontre que la vaccination en anneau constitue une stratégie d'intervention viable et hautement efficiente pour contrôler la propagation d'un pathogène mortel sous des contraintes budgétaires strictes. En mobilisant un ciblage spatial précis dès la détection du premier décès, cette approche permet de réduire les coûts de 77,3 % par rapport à une vaccination de masse, tout en assurant un taux de survie moyen de 45,1 % (soit 1 692 individus sauvés) face à une maladie 100 % létale.

L'analyse des 100 réplications souligne toutefois que le succès de cette stratégie est étroitement lié à la stochasticité spatiale. Comme l'illustre la distribution bimodale des résultats, l'efficacité de l'intervention dépend de la rapidité avec laquelle les premières barrières immunitaires parviennent à isoler les foyers infectieux. Ce modèle confirme qu'une gestion optimisée des ressources (tests RAT et doses de vaccins) peut pallier l'absence d'information parfaite sur la prévalence épidémique, à condition que la riposte soit géographiquement localisée.

## Limites du modèle

Malgré la pertinence des résultats obtenus, plusieurs simplifications limitent l'application directe de ce modèle à des contextes réels :

- **Dynamique biologique** : L'absence d'une période d'incubation masque la phase de transmission invisible. De plus, l'hypothèse d'une immunité vaccinale absolue et permanente ne reflète pas la réalité de la réponse immunitaire qui nécessiterait des doses de rappel (booster shots).
- **Hétérogénéité spatiale et comportementale** : La population est répartie uniformément sur la grille. L'intégration de zones de haute densité (hubs comme les écoles ou bureaux) et de comportements adaptatifs permettrait de mieux capturer les phénomènes de super-propagation.
- **Contraintes administratives** : L'impossibilité d'intervenir avant le premier décès crée une latence critique. Dans un scénario réel, la variabilité de la durée de la maladie et de la létalité influencerait grandement le timing de la riposte.
