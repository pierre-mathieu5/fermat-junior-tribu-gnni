# Reconstruire une tribu : une géométrie cachée

Programme qui accompagne le manuscrit.

Il refait les calculs de la Section 7 du mémoire : pour n = 3, …, 7, le nombre de sommets et d'arêtes du graphe H_n et son diamètre (Tableau 1), ainsi que le décompte des exceptions pour n = 3 et n = 4.

## Utilisation

Il faut Python 3.8 ou plus récent, avec `numpy` et `scipy` :

```
pip install numpy scipy
python3 calculs_tribu.py
```

Le calcul dure environ une minute et affiche :

```
n = 3 : 16 sommets, 48 aretes, diametre 3 (4 formes)
        0 exceptions, dont 0 provenant d'une transposition
n = 4 : 125 sommets, 730 aretes, diametre 4 (9 formes)
        30 exceptions, dont 30 provenant d'une transposition
n = 5 : 1296 sommets, 12120 aretes, diametre 5 (20 formes)
n = 6 : 16807 sommets, 225141 aretes, diametre 7 (48 formes)
n = 7 : 262144 sommets, 4682496 aretes, diametre 9 (115 formes)
```

## Méthode

Les arbres étiquetés sont énumérés par leurs codes de Prüfer, puis enracinés en \* pour obtenir leurs familles de clusters. Deux arbres sont voisins dans H_n lorsque leurs familles ont n − 1 clusters en commun. Les distances sont calculées par des parcours en largeur.

Pour le diamètre, nous n'avons pas besoin de partir de chaque arbre. Renommer les points de B envoie les arêtes de H_n sur des arêtes, donc l'excentricité d'un arbre ne dépend que de sa forme, c'est-à-dire de l'arbre enraciné non étiqueté sous-jacent. Un parcours en largeur depuis un arbre de chaque forme suffit : il y en a 115 pour n = 7.
