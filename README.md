# Electoral Space

Para representar 3 variables en 2 dimensiones de forma única para cada conjunto de valores proporcionales se puede utilizar un diagrama ternario. Se puede entender el diagrama como la proyección de los puntos del cuadrante positivo del espacio en el hiperplano cuya suma de coordenadas es 1. La imagen de la proyección es un triángulo equilátero. De la misma manera podemos representar dos valores en un segmento o cuatro valores en un tetraedro. Esta representación también se conoce con el nombre de coordenadas baricéntricas en un símplice.

Un diagrama como éste se puede usar para representar el resultado de unas elecciones en las que los votos se reparten entre 3 partidos. Cada punto del triángulo se representa con tres coordenadas que se corresponden con el porcentaje de voto obtenido por cada partido. Tenemos un número de escaños a repartir y cada posible reparto se representa con un punto resaltado (nodo) en el triángulo con su correspondiente resultado. El método electoral asignará a cada resultado electoral uno de los posibles repartos de escaños. Las regiones de puntos a las que le son asignados un mismo reparto serán representadas con un mismo color. Para dos partidos usaríamos un segmento y para cuatro, un tetraedro. El diagrama es útil para visualizar los posibles escaños en juego cuando el resultado está muy ajustado y próximo a la frontera entre dos o más regiones.

Utilizamos una función de asignación electoral para distintos métodos proporcionales de cocientes con la posibilidad de añadir un umbral de entrada. 

Mostramos algunas de las funciones utilizadas para calcular la asignación y generar los diagramas.

Las regiones creadas por el reparto de D'Hondt no se corresponden con las regiones de Voronoi de los nodos. Este método asigna correctamente más del 60% de los puntos. El reparto de Sainte-Laguë se aproxima en más de un 90% a las regiones de Voronoi sin llegar a ser igual. Con D'Hondt el tamaño de las distintas regiones es similar, mientras para Sainte-Laguë los nodos se encuentran en el centro de ellas. Obtendríamos las regiones de Voronoi si a cada punto se le asignara el nodo más próximo. Para este tipo de reparto se obtienen los mismos resultados si utilizamos la distancia Euclídea, Manhattan o la uniforme. Como detalle, el cálculo para la asignación de Voronoi debería resolver un empate entre dos o más partidos cuando el resultado electoral se situa exactamente en la frontera entre dos o más regiones contiguas, del mismo modo que los métodos de cocientes resuelve los empates de cocientes (con el número de votos total y al azar en caso de nuevo empate).

Los métodos de cocientes producen regiones con fronteras dadas por sectores angulares mientras que la asignación de Voronoi produce regiones dadas por sectores de segmentos.

Se estudia aquí la asignación no sólo en resultado, sino también en el orden de reparto de los escaños, y se dibujan diagramas que muestran las regiones del Espacio Electoral divididas en subregiones para cada posible ordenación de la asignación. El diagrama de las subregiones visibiliza de manera bastante destacable la geometría de las regiones electorales.

Podemos también representar datos históricos o geográficos de varias elecciones sobre el diagrama.

En http://rpubs.com/ismaelsb/ElectoralSpace o en https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.md se presentan los resultados 

El código utilizado para estos cálculos y representaciones se encuentra en este enlace https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.R

Los resultados han sido generados en HTML a través del código R Markdown https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.Rmd
