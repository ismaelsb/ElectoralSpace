﻿# ElectoralSpace

Este diagrama ternario representa el resultado de unas elecciones en las que los votos se reparten entre 3 partidos. Cada punto del triángulo se representa con tres coordenadas que se corresponden con el porcentaje de voto obtenido por cada partido. Tenemos un número de escaños a repartir y cada posible reparto se representa con un punto negro (nodo) en el triángulo con su correspondiente resultado. El método electoral asignará a cada resultado electoral uno de los posibles repartos de escaños. Las regiones de puntos a las que le son asignados un mismo reparto están representadas con un mismo color. Las regiones creadas por el reparto de D'Hondt no se corresponden con regiones de Voronoi de los nodos. Asigna correctamente más del 60% de los puntos. El reparto de Sainte-Laguë se aproxima más, alrededor de un 90%, a las regiones de Voronoi pero no llega a ser igual. En D'Hondt el tamaño de las distintas regiones es similar, en Sainte-Laguë los nodos se encuentran en el centro de ellas. Las regiones de Voronoi resultarían si a cada punto se le asignara el nodo más próximo. Para este tipo de reparto se obtienen los mismos resultados si utilizamos la distancia Euclídea, Manhattan o la uniforme. El triángulo representa la región del espacio ΣXi=1, 0≤Xi≤1. Para 2 partidos usaríamos un segmento y para 4, un tetraedro. El diagrama es útil para visualizar los posibles escaños en juego cuando el resultado está muy ajustado y próximo a la frontera entre dos o más regiones.
El código utilizado para estos cálculos y representaciones se encuentra en este enlace https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.R

En el código se incluyen la función de asignación electoral y la asignación con umbral de entrada y se muestran los efectos del umbral electoral sobre las regiones de asignación en un diagrama. Se añade además la representación de datos históricos o geográficos de varias elecciones sobre el diagrama. Se estudia la asignación no sólo en resultado, sino también en el orden de reparto de los escaños, y se dibujan diagramas que muestran las regiones del Espacio Electoral divididas en subregiones para cada posible ordenación de la asignación. Uno de los diagramas de las subregiones visibiliza de manera bastante destacable la geometría de las regiones electorales.

Se puede ver una muestra de los resultados en https://dl.dropboxusercontent.com/u/1570957/ElectoralSpace.html que han sido generados en HTML a través del código R Markdown https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.Rmd
