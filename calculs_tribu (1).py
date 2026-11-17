import itertools
from collections import defaultdict

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import shortest_path

# Les points de B sont numérotés 0, ..., n-1, et la racine * porte le numéro n.


def arbres(n):
    """Arbres étiquetés sur B u {*}, énumérés par leurs codes de Prüfer."""
    N = n + 1
    for code in itertools.product(range(N), repeat=N - 2):
        # Nous décodons la suite de Prüfer : à chaque étape, la plus petite feuille
        # est reliée au terme courant de la suite, puis retirée.
        deg = [1] * N
        for a in code:
            deg[a] += 1
        aretes = []
        for a in code:
            f = min(i for i in range(N) if deg[i] == 1)
            aretes.append((f, a))
            deg[f] -= 1
            deg[a] -= 1
        u, v = [i for i in range(N) if deg[i] == 1]
        aretes.append((u, v))
        yield aretes


def clusters(n, aretes):
    """Famille de clusters de l'arbre enraciné en *, et sa forme (arbre non étiqueté)."""
    voisins = defaultdict(list)
    for a, b in aretes:
        voisins[a].append(b)
        voisins[b].append(a)
    # Parcours depuis * pour connaître le parent de chaque sommet.
    parent = {n: None}
    ordre = [n]
    for v in ordre:
        for w in voisins[v]:
            if w not in parent:
                parent[w] = v
                ordre.append(w)
    enfants = defaultdict(list)
    for v in ordre[1:]:
        enfants[parent[v]].append(v)
    # Le cluster S_v est formé de v et de ses descendants ; nous le codons par un entier
    # dont le bit i vaut 1 si i est dans S_v. On remonte depuis les feuilles.
    S, forme = {}, {}
    for v in reversed(ordre):
        S[v] = (1 << v) | sum(S[c] for c in enfants[v])
        forme[v] = "(" + "".join(sorted(forme[c] for c in enfants[v])) + ")"
    return frozenset(S[v] for v in range(n)), forme[n]


def graphe(n):
    """Le graphe H_n : deux familles sont voisines si elles ont n-1 clusters en commun."""
    familles, formes = [], []
    for aretes in arbres(n):
        F, forme = clusters(n, aretes)
        familles.append(F)
        formes.append(forme)
    # Deux familles voisines ont la même famille une fois retiré le cluster qui change :
    # il suffit de regrouper les familles selon chaque sous-famille à n-1 clusters.
    groupes = defaultdict(list)
    for i, F in enumerate(familles):
        for C in F:
            groupes[F - {C}].append(i)
    I, J = [], []
    for g in groupes.values():
        for a, b in itertools.combinations(g, 2):
            I.append(a)
            J.append(b)
    M = len(familles)
    A = coo_matrix((np.ones(2 * len(I), dtype=np.int8), (I + J, J + I)), shape=(M, M)).tocsr()
    return familles, formes, A, len(I)


def echange(F, x, y):
    """Clusters de l'arbre obtenu en échangeant les étiquettes x et y."""
    G = set()
    for C in F:
        bx, by = (C >> x) & 1, (C >> y) & 1
        G.add(C & ~((1 << x) | (1 << y)) | (bx << y) | (by << x))
    return frozenset(G)


for n in range(3, 8):
    familles, formes, A, nb_aretes = graphe(n)
    assert len(familles) == (n + 1) ** (n - 1)  # formule de Cayley

    # Renommer les points de B ne change pas les distances : l'excentricité d'un arbre
    # ne dépend donc que de sa forme, et il suffit de partir d'un arbre de chaque forme.
    representants = list({f: i for i, f in reversed(list(enumerate(formes)))}.values())
    D = shortest_path(A, unweighted=True, directed=False, indices=representants)
    assert np.isfinite(D).all()  # le graphe est connexe
    print(
        f"n = {n} : {len(familles)} sommets, {nb_aretes} aretes, "
        f"diametre {int(D.max())} ({len(representants)} formes)"
    )

    # Exceptions : paires de sommets pour lesquelles la borne triviale n'est pas atteinte.
    if n <= 4:
        Dtout = shortest_path(A, unweighted=True, directed=False)
        indice = {F: i for i, F in enumerate(familles)}
        transpositions = set()
        for i, F in enumerate(familles):
            for x, y in itertools.combinations(range(n), 2):
                j = indice[echange(F, x, y)]
                transpositions.add((min(i, j), max(i, j)))
        exceptions = [
            (i, j)
            for i, j in itertools.combinations(range(len(familles)), 2)
            if Dtout[i, j] > n - len(familles[i] & familles[j])
        ]
        dans = sum(1 for e in exceptions if e in transpositions)
        print(f"        {len(exceptions)} exceptions, dont {dans} provenant d'une transposition")
