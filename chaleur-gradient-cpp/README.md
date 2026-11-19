# Équation de la chaleur 1D : méthodes de gradient en C++

Projet d'Analyse Matricielle (M1 Mathématiques et Interactions, Paris-Saclay Évry).

On résout l'équation de la chaleur 1D par un schéma d'Euler implicite. À chaque pas de temps, on doit résoudre un système linéaire `A x = b` où `A` est symétrique définie positive (tridiagonale). Ce système est résolu par trois méthodes itératives :

- gradient à **pas fixe**
- gradient à **pas optimal**
- gradient **conjugué**

La matrice peut être stockée de deux façons : **stockage creux** (type COO : lignes, colonnes, valeurs) ou **stockage classique** (matrice dense). L'erreur max est calculée par rapport à la solution exacte `u(t,x) = (1+16384 µ t)^(-1/2) exp(-4096 (x-1/2)² / (1+16384 µ t))`.

## Structure

```
chaleur-gradient-cpp/
├── include/
│   ├── ressource.hpp   # classe Stockage (COO), outils vectoriels, erreur
│   └── gradient.hpp    # gradientPasFixe, gradient_optimal, gradientConjugue
├── src/
│   ├── main.cpp        # programme principal (ancien test.cpp)
│   ├── fonction.cpp    # matrice, produits, fonctionnelle
│   └── erreur.cpp      # condition initiale, solution exacte, erreur max
├── docs/
│   └── rendu_ecrit.pdf # partie théorique + tests et interprétations
├── results/
│   ├── resultats_exemple.txt   # sortie (x, u(x)) d'une exécution
│   └── figures/                # graphiques par méthode et type de stockage
├── Makefile
└── README.md
```

## Compilation et exécution

```bash
make
./chaleur
```

Le programme demande le type de stockage (0 creux, 1 classique) puis la méthode (0 pas fixe, 1 optimal, 2 conjugué). Il écrit `resultats.txt` (positions et valeurs de la solution) et affiche le temps final et l'erreur max.

Les paramètres (`N`, `alpha`, `Tf`, `iter`) se règlent en haut de `src/main.cpp`. Pour les tests du rapport : `Tf = 0.02`, `iter = 1000`, `alpha = 0.5`.

## Résultats (résumé du rapport)

Pour `Tf = 0.02`, 1000 itérations, `alpha = 0.5` : seul le gradient conjugué atteint une erreur max de l'ordre de 1e-5. Le meilleur compromis précision / temps est `(alpha, N) = (0.5, 256)`. Voir `docs/rendu_ecrit.pdf`.

## Auteur

Clovis Retif
