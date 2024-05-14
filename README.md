# Éléments Finis en C++
- main.cpp :
  Programme principal permettant de lancer les tests ainsi que les simulations.
  
- src/tests.h :
Contient la définition des différents tests qui seront appelés dans main.cpp. 
- Les différentes tests sont :
  test-quadrature : permet de tester la structure Quadrature
  test_elemnt_mapping : teste les méthodes de la classe ElementMapping, qui permet de transfomer les coordonnées de l'élement de réference vers l'espace global, calculer le jacobien du mapping ainsi que son déterminant
  test_shape_functions : teste la classe ShapeFunctions, et renseigne le nombre de fonctions de forme en fonction de l'élément, l'expression des fonctions de forme ainsi que leur gradient
  test_assemble_elementary_matrix : teste la méthode assemble_elementary_matrix, qui construit Ke
  test_assemble_elementary_vector : teste la méthode assemble_elementary_vector, qui construit Fe

- src/simu.h :
Contient la définition des différentes simulations qui seront appelées dans
main.cpp.

- src/fem.cpp :
Contient l'implémentation des méthodes des éléments finis.

  
#### Contact

Paul Cupillard : paul.cupillard@univ-lorraine.fr

#### Liens utiles

Cours sur [Arche](http://arche.univ-lorraine.fr/course/view.php?id=61482)

Vidéo de Gilbert Strang : [Finite element method](https://www.youtube.com/watch?v=WwgrAH-IMOk)

Cours de Grégoire Allaire : [Approximation numérique et optimisation](http://www.cmap.polytechnique.fr/~allaire/map411/polycopie-map411.pdf)

Générateurs de maillages triangulaires : [Gmsh](http://gmsh.info/),[Triangle](https://www.cs.cmu.edu/~quake/triangle.html)
