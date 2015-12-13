# Electoral Space. Geometric Visualizacion of Electoral Proportional Methods.


In order to represent 3 variables in 2 dimensions in an only way for each set of proportional values we can use a ternary diagram. You can see this diagram as the projection of the points from the positive quadrant in space onto the hyperplane of points whose sum of coodinates equals 1. The image of this projection is an equilateral triangle. In the same way we can represent two values in a segment and four values in a tetrahedron. This representation is also known by the name of barycentric coordinates in a simplex.
 
A diagram like these can be used to represent the results in political elections in which the votes are shared between 3 parties. Each point in the triangle is represented by three coordinates corresponding the share of votes received by each party. We have a determined number of seats to allocate, and each posible sharing is represented by a highlighted dot (a node) in the triangle and a label with its result. The electoral method would allocate to each electoral result one of the possible sharings of seats. The regions of points to which the same sharing is allocated are represented in the same color. For two parties we'd use a segment. For four, a tetrahedron. The diagram is usefull for visualizing possible seats at stake when the results are close to the borders between two or more regions.

We use a function for electoral allocation for different proportional methods and an optional entry threshold.

We also show some of the functions used to compute allocations or to plot the diagrams.

The regions created by D'Hondt method aren't the same as the Voronoi regions created by the lattice of nodes. This method allocates more than 60% of the points correctly. Sainte-Laguë method is, in more than 90% of the points, similar to Voronoi. In D'Hondt the sizes of the regions are similar, but Sainte-Laguë produces regions centered on the nodes. Voronoi regions would be obtained by a method always allocating the nearest node. Those regions are similar if we use Euclidean, Manhattan or uniform distance. A Voronoi allocation method should solve the ties between two or more parties when the result is in a border in the same way the quotients methods do (by the total number of votes, and if the tie persists, at random). We study two types of Voronoi regions: those created by distances on the projective space of the simplex (Euclidean, Manhattan and Uniform give the same regions), and those created by the orthodromic distance on the unit sphere (normalized vectors).

Quotients method produce regions with borders in angular sections. Voronoi borders are given by segment sections.

We study here not only the result of allocation, but also the ordering in the allocation, and so we plot a diagram for the regions of different allocation orderings in the Electoral Space as subregions of the allocation ones.
This diagram visualizes in a remarkable way the geometry of the Electoral Space.

We can also draw historic data from past elections on the diagrams.

Next we show vote disperion and camera dispersion measured with different diversity indexes: Shannon entropy and Laakso-Taagepera effective number of parties.

It can be observed that the Largest remainder method with Hare Quota gives the nearest allocation to each share of votes. Its allocation regions are the Voronoi regions given by the seat allocation nodes.

Hare quota is pretty close to Sainte-Laguë and Droop quota is so to D'Hondt as can be seen in the next diagrams.

In http://rpubs.com/ismaelsb/ElectoralSpace and https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.md we present the results.

The code used for these computations and representations is available here: https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.R

The results have been generated in HTML by this R Markdown code: https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.Rmd

An interactive web app is available here: https://ismaelsb.shinyapps.io/ElectoralSpace





Para representar 3 variables en 2 dimensiones de forma única para cada conjunto de valores proporcionales se puede utilizar un diagrama ternario. Se puede entender el diagrama como la proyección de los puntos del cuadrante positivo del espacio en el hiperplano cuya suma de coordenadas es 1. La imagen de la proyección es un triángulo equilátero. De la misma manera podemos representar dos valores en un segmento o cuatro valores en un tetraedro. Esta representación también se conoce con el nombre de coordenadas baricéntricas en un símplice.

Un diagrama como éste se puede usar para representar el resultado de unas elecciones en las que los votos se reparten entre 3 partidos. Cada punto del triángulo se representa con tres coordenadas que se corresponden con el porcentaje de voto obtenido por cada partido. Tenemos un número de escaños a repartir y cada posible reparto se representa con un punto resaltado (nodo) en el triángulo con su correspondiente resultado. El método electoral asignará a cada resultado electoral uno de los posibles repartos de escaños. Las regiones de puntos a las que le son asignados un mismo reparto serán representadas con un mismo color. Para dos partidos usaríamos un segmento y para cuatro, un tetraedro. El diagrama es útil para visualizar los posibles escaños en juego cuando el resultado está muy ajustado y próximo a la frontera entre dos o más regiones.

Utilizamos una función de asignación electoral para distintos métodos proporcionales de cocientes con la posibilidad de añadir un umbral de entrada. 

Mostramos algunas de las funciones utilizadas para calcular la asignación y generar los diagramas.

Las regiones creadas por el reparto de D'Hondt no se corresponden con las regiones de Voronoi de los nodos. Este método asigna correctamente más del 60% de los puntos. El reparto de Sainte-Laguë se aproxima en más de un 90% a las regiones de Voronoi sin llegar a ser igual. Con D'Hondt el tamaño de las distintas regiones es similar, mientras para Sainte-Laguë los nodos se encuentran en el centro de ellas. Obtendríamos las regiones de Voronoi si a cada punto se le asignara el nodo más próximo. Para este tipo de reparto se obtienen los mismos resultados si utilizamos la distancia Euclídea, Manhattan o la uniforme. Como detalle, el cálculo para la asignación de Voronoi debería resolver un empate entre dos o más partidos cuando el resultado electoral se situa exactamente en la frontera entre dos o más regiones contiguas, del mismo modo que los métodos de cocientes resuelve los empates de cocientes (con el número de votos total y al azar en caso de nuevo empate). Estudiamos dos tipos de regiones de Voronoi: las creadas por una distancia en el espacio proyectivo del símplice (la euclídea, la Manhattan y la uniforme dan las mismas regiones), y las creadas por la distancia ortodrómica en la esfera unidad (con vectores normalizados).

Los métodos de cocientes producen regiones con fronteras dadas por sectores angulares mientras que la asignación de Voronoi produce regiones dadas por sectores de segmentos.

Se estudia aquí la asignación no sólo en resultado, sino también en el orden de reparto de los escaños, y se dibujan diagramas que muestran las regiones del Espacio Electoral divididas en subregiones para cada posible ordenación de la asignación. El diagrama de las subregiones visibiliza de manera bastante destacable la geometría de las regiones electorales.

Podemos también representar datos históricos o geográficos de varias elecciones sobre el diagrama.

A continuación mostramos la dispersión de los votos y de la cámara medidas con diferentes índices de diversidad: la entropía de Shannon y el número efectivo de partidos de Laakso-Taagepera.

Se puede observar que el método de los restos mayores con el cociente Hare da la asignación más cercana a cada reparto de votos. Su regiones de asignación son las regiones de Voronoi dadas por los nodos de asignación de escaños.

El cociente Hare es bastante parecido a Sainte-Laguë y el cociente Droop lo es a D'Hondt, como puede verse en los siguientes diagramas.

En http://rpubs.com/ismaelsb/ElectoralSpace o en https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.md se presentan los resultados 

El código utilizado para estos cálculos y representaciones se encuentra en este enlace https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.R

Los resultados han sido generados en HTML a través del código R Markdown https://github.com/ismaelsb/ElectoralSpace/blob/master/ElectoralSpace.Rmd
