# Electoral Space

Código disponible en <https://github.com/ismaelsb/ElectoralSpace>
 


**Introducción**

Un diagrama ternario se utiliza para representar 3 variables en 2 dimensiones de forma única para cada conjunto de valores proporcionales. Se puede entender el diagrama como la proyección de los puntos del cuadrante positivo del espacio en
el hiperplano cuya suma de coordenadas es 1. La imagen de la proyección es un triángulo equilátero. De la misma manera podemos representar dos valores en un segmento o cuatro valores en un tetraedro.



![](ElectoralSpace_files/figure-html/unnamed-chunk-2-1.png) 


**Funciones de asignación electoral (más en el código)**

Un diagrama como éste se puede usar para representar el resultado de unas elecciones en las que los votos se reparten entre 3 partidos. Cada punto del triángulo se representa con tres coordenadas que se corresponden con el porcentaje de voto obtenido por cada partido. Tenemos un número de escaños a repartir y cada posible reparto se representa con un punto resaltado (nodo) en el triángulo con su correspondiente resultado. El método electoral asignará a cada resultado electoral uno de los posibles repartos de escaños. Las regiones de puntos a las que le son asignados un mismo reparto están representadas con un mismo color. Para dos partidos usaríamos un segmento y para cuatro, un tetraedro. El diagrama es útil para visualizar los posibles escaños en juego cuando el resultado está muy ajustado y próximo a la frontera entre dos o más regiones.

Utilizamos una función de asignación electoral para distintos métodos proporcionales de cocientes con la posibilidad de añadir un umbral de entrada. 

Mostramos algunas de las funciones utilizadas para calcular la asignación y generar los diagramas:


```r
alloc <- function(parties, votes, seats, step, threshold=0){
  
  #function for the seat allocation and its ordering
  #computes allocation for several values of seats and divisor step methods
  #with the same quotiens table
  votes=votes*(votes>=(threshold*sum(votes)))
  Mst <- max(step)
  Mse <- max(seats)
  nquotients=(1+Mst*(Mse-1))
  
  #table with all the quotients needed
  quotienstable <- data.frame( 
    parties    = rep(parties, each = nquotients), 
    quotients  = as.vector(sapply(votes, function(x) 
      x/seq(from=1, to=nquotients) )),
    votesrep   = rep(votes, each = nquotients)
  ) 
  
  SeatsList=list()
  filteredtable =list()
  
  
  for (j in 1:length(step)) {
    
    
    #select from the table of all quotiens the ones whose divisors
    #belong to the current sequence, given 'seats' and 'step'
    select <- rep(seq(from=1, to=1+step[j]*(Mse-1), by=step[j]),length(parties))+
      rep(seq(from=0,to=nquotients*length(parties)-1,by=nquotients),each=Mse)
    
    filteredtable [[j]] <- quotienstable$parties[select][order(-quotienstable$quotients[select], -quotienstable$votesrep[select])]
    
    SeatsList[[j]]<-list()
    
    #a vector of values for seats if you want to compute the partial sums
    for (i in 1:length(seats)){
      SeatsList[[j]][[i]]= table(filteredtable[[j]][1:seats[i]])
    }
    
    names(SeatsList[[j]])[1:length(seats)]<-do.call("paste0",as.data.frame(cbind("divisor step ",step[j]," for ", seats," seats")))
    
    
    SeatsList[[j]][[length(seats)+1]] = as.matrix(filteredtable[[j]][1:Mse]) #ordering
    
    names(SeatsList[[j]])[length(seats)+1]<-paste0("ordering for divisor step ", step[j]," for ", Mse," seats")
    
    
  }
  
  
  return(SeatsList);
  
}
```


```r
generateNodes <- function(seats){ 
  
  #define nodes
  nnodes = (seats+1)*(seats+2)/2;
  nodes  = matrix(0,nnodes,3);
  t=1;
  for (i in 0:seats){
    for (j in i:seats){
      nodes[t,] = c(i,j-i,seats-j);
      t=t+1;
    }
  }   
  return(nodes);
  
}

generateDots <- function(dots, method="cartesian"){
  
  if (method == "cartesian") {
    #Map Cartesian to Ternary to produce a homegeneous simulation
    Rc = matrix(runif(2*dots), nrow=dots, ncol=2)
    R  = matrix(apply(Rc,1, CartesianToTernary), nrow=dots, ncol=3, byrow=TRUE)
  }
  
  else if (method == "ternary") {
    #Ternary simulation. Projects points in[0,1]^3 to the sum(x)=1 hyperplane
    #This creates way more points in the center than in the extremes
    R = matrix(runif(3*dots), nrow=dots, ncol=3)
    R = prop.table(R,1) #rows sum 1
  }
  
  return(R);
  
}
```




```r
ManhattanNearest <- function (x, seats, nodes, nnodes) {
  
  Manhattan = which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,sum))
  
  return(Manhattan);
  
}
```




```r
generateSpline <- function (dfvotes, method = "natural") {
  
  #generates a spline curve through elections history
  
  nelect = dim(dfvotes)[1]
  
  splineX <- spline(x=1:nelect,y=dfvotes$x, method=method)
  splineY <- spline(x=1:nelect,y=dfvotes$y, method=method)
  
  dfSpline = matrix(0, length(splineX[[1]]), 3)
  dfSpline[,1] <- splineX[[2]]
  dfSpline[,2] <- splineY[[2]]
  dfSpline[,3] <- 1-dfSpline[,1]-dfSpline[,2]
  dfSpline = as.data.frame(dfSpline)
  names(dfSpline) <- letters[24:26]
  
  return(dfSpline);
  
}

generateColors <- function (colorRGB, seats) {
  
  nnodes = (seats+1)*(seats+2)/2;
  nodes <- generateNodes(seats)
  #decimal codes for colors in nodes
  Code3 = floor(nodes*255.9/seats) #255.9 avoids seats -> 256 -> HEX #100 case
  #linear transformation to fixed extreme colors
  CodeDecRGB = floor(Code3%*%matrix(colorRGB,3,3, byrow=T)/256)
  #hex codes for colors in nodes
  CodeRGB=cbind(matrix("#",nnodes),format(as.hexmode(CodeDecRGB[,1]),width=2),format(as.hexmode(CodeDecRGB[,2]),width=2),format(as.hexmode(CodeDecRGB[,3]),width=2))
  CodeRGB=do.call("paste0",as.data.frame(CodeRGB))
  values=CodeRGB
  
  return(values);
  
}

generateOrderColors <- function (colorRGB, seats, validCodes=TRUE) {

  #generate order color from alloc order codes
  
  #extract allocation order from allocation order codes:
  OrderCode <-as.matrix(0:(3^seats-1))

  OrderRests <- t(matrix(apply(OrderCode, 1, function(x) x %% (3^(1:seats))),ncol=3^seats, nrow=seats))
  #(this is vectorial for OrderRests <- OrderCode %% (3^(1:seats))  )
  
  Order_from_code <- (cbind(OrderRests,0)-cbind(0,OrderRests))[,1:seats]/t(matrix(rep(3^(0:(seats-1)),3^seats),ncol=3^seats,nrow=seats))+1
  vectorOrderCode3 <- 2^(0:(seats-1))
  #decimal codes for colors in nodes
  Code3 = cbind((Order_from_code == 1)%*%as.matrix(vectorOrderCode3),
                (Order_from_code == 2)%*%as.matrix(vectorOrderCode3),
                (Order_from_code == 3)%*%as.matrix(vectorOrderCode3))
  
  #linear transformation to fixed extreme colors
  CodeDecRGB = floor(Code3%*%matrix(colorRGB,3,3, byrow=T)/2^seats)
  
  #hex codes for colors in nodes
  CodeRGB = cbind(matrix("#",3^seats),format(as.hexmode(CodeDecRGB[,1]),width=2),format(as.hexmode(CodeDecRGB[,2]),width=2),format(as.hexmode(CodeDecRGB[,3]),width=2))
  CodeRGB = do.call("paste0",as.data.frame(CodeRGB))

  #values = CodeRGB[validCodes]
  values = CodeRGB[c(validCodes, 3^seats)] #added one more color bc scale_color_manual() demanded one more color
  #values = CodeRGB[c(1, validCodes)] #added one more color bc scale_color_manual() demanded one more color
  
  
return(values);

}
```


**Matriz de datos (más en el código)**


```r
AssembleData <- function (dots, seats, step=1, threshold=0) {
  
  
  #Generate random electoral results
  R <- generateDots(dots)
  
  Seats=max(seats)
  
  #nodes
  nodes  = generateNodes(Seats)
  nnodes = (Seats+1)*(Seats+2)/2;
  
  #node labels for seats
  label=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
  label=do.call("paste0",as.data.frame(label))
  label=rbind(matrix("",dots),as.matrix(label))
  
  #attach node points                
  x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/Seats))
  y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/Seats))
  z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/Seats))
  
  #Indexes for Voronoi regions
  Uniform   = apply(R, 1, UniformNearest,   seats=Seats, nodes = nodes, nnodes=nnodes)
  Manhattan = apply(R, 1, ManhattanNearest, seats=Seats, nodes = nodes, nnodes=nnodes)
  Euclid    = apply(R, 1, EuclidNearest,    seats=Seats, nodes = nodes, nnodes=nnodes)
  
  #allocate seats
  AllocStructure <- apply(R, 1, function(x) alloc(as.character(1:3), x, seats, step, threshold))
  #input a vector of values for seats of to compute partial sums
  
  df=list() #list of dataframes, one for each 'step' value
  
  #loop for diferent 'step's
  for (j in 1:length(step)){
    
    #allocation
    AllocPartial=matrix(0,dots,length(seats))
    
    for (i in 1:length(seats)){
      
      S = t(matrix(sapply(AllocStructure, function(x) x[[j]][[i]]),nrow=3,ncol=dots))
      
      AllocPartial[,i] = apply(S, 1, AllocatedNode, nodes = generateNodes(seats[i]), nnodes=(seats[i]+1)*(seats[i]+2)/2)
      
    }
    
    #allocation order
    #computes allocation ordering for max(seats) in each step value
    AllocOrder = t(matrix(sapply(AllocStructure, function(x) as.integer(x[[j]][[length(seats)+1]])),nrow=Seats,ncol=dots))
    vectorOrderCode=as.matrix(3^(0:(Seats-1)))
    #this coding highlights the last seats over the first ones
    # so that there is contrast between adjacent regions
    AllocOrderCode = (matrix(AllocOrder, ncol=Seats)-1) %*% as.matrix(vectorOrderCode)
    
    AllocPartial=rbind(AllocPartial,matrix(NA,nnodes,length(seats)))
    dfPartial=as.data.frame(AllocPartial)
    names(dfPartial) <- do.call(paste0,as.data.frame(cbind("All",seats)))
    names(dfPartial)[which(seats==Seats)] <- "Allocated"
    
    #assembling a data frame
    df0 = data.frame(
      
      type          = c(matrix("dot",dots), matrix("node",nnodes)),
      
      x ,
      y ,
      z ,
      
      Sx            = c(S[,1], matrix(NA,nnodes)),
      Sy            = c(S[,2], matrix(NA,nnodes)),
      Sz            = c(S[,3], matrix(NA,nnodes)),
      
      Euclid        = c(Euclid, matrix(NA,nnodes)),
      Manhattan     = c(Manhattan, matrix(NA,nnodes)),
      Uniform       = c(Uniform, matrix(NA,nnodes)),
      
      Malapportionment = c(AllocPartial[1:dots,length(seats)] != Euclid[1:dots], matrix(NA,nnodes)),
      
      AllocOrderCode=  c(AllocOrderCode, matrix(NA,nnodes)),
      
      label
      
    )
    
    df0=cbind(df0,dfPartial)
    
    df[[j]]<-df0
    
  } #end of loop for different 'step's
  
  
  return(df);
  
}
```




**Ejemplo de asignación**


```
## [1]  77 265 117
```

```
## [[1]]
## [[1]]$`divisor step 1 for 9 seats`
## 
## a b c 
## 1 6 2 
## 
## [[1]]$`ordering for divisor step 1 for 9 seats`
##       [,1]
##  [1,] "b" 
##  [2,] "b" 
##  [3,] "c" 
##  [4,] "b" 
##  [5,] "a" 
##  [6,] "b" 
##  [7,] "c" 
##  [8,] "b" 
##  [9,] "b" 
## 
## 
## [[2]]
## [[2]]$`divisor step 2 for 9 seats`
## 
## a b c 
## 2 5 2 
## 
## [[2]]$`ordering for divisor step 2 for 9 seats`
##       [,1]
##  [1,] "b" 
##  [2,] "c" 
##  [3,] "b" 
##  [4,] "a" 
##  [5,] "b" 
##  [6,] "c" 
##  [7,] "b" 
##  [8,] "b" 
##  [9,] "a"
```


**Configuración**


```r
#presets
seats=2:5;
step=c(1,2); #(2 Sainte-Laguë 1 D'Hondt)
dots=40000
threshold=0
```


**Generando los datos**


```r
#Assemble the data frame

df = AssembleData(dots, seats, step)

#df = AssembleData(dots, 5, 1, threshold)

#df = AssembleData(dots, c(3,5,4), c(2,1), threshold)

dfT = AssembleData(dots, seats=5, threshold=.20)
```




**Regiones del Espacio Electoral**

Las regiones creadas por el reparto de D'Hondt no se corresponden con regiones de Voronoi de los nodos. Asigna correctamente más del 60% de los puntos. El reparto de Sainte-Laguë se aproxima más, más de un 90%, a las regiones de Voronoi pero no llega a ser igual. En D'Hondt el tamaño de las distintas regiones es similar, en Sainte-Laguë los nodos se encuentran en el centro de ellas. Las regiones de Voronoi resultarían si a cada punto se le asignara el nodo más próximo. Para este tipo de reparto se obtienen los mismos resultados si utilizamos la distancia Euclídea, Manhattan o la uniforme.

Se estudia aquí la asignación no sólo en resultado, sino también en el orden de reparto de los escaños, y se dibujan diagramas que muestran las regiones del Espacio Electoral divididas en subregiones para cada posible ordenación de la asignación. El diagrama de las subregiones visibiliza de manera bastante destacable la geometría de las regiones electorales.

![](ElectoralSpace_files/figure-html/unnamed-chunk-15-1.png) ![](ElectoralSpace_files/figure-html/unnamed-chunk-15-2.png) 


**Comparación de las regiones de Voronoi sobre distintas métricas**

Las distancias euclídea, Manhattan y uniforme producen regiones de Voronoi idénticas:


```
## [1] 1
```

```
## [1] 1
```

```
## [1] 1
```


**Tamaño de las regiones**

![](ElectoralSpace_files/figure-html/unnamed-chunk-17-1.png) ![](ElectoralSpace_files/figure-html/unnamed-chunk-17-2.png) 


**Proporción de resultados no asignados al nodo más cercano**


```
## [1] 0.365225
```

```
## [1] 0.075875
```


**Efecto del umbral sobre las regiones**

![](ElectoralSpace_files/figure-html/unnamed-chunk-19-1.png) 


**Regiones para las distintas ordenaciones en la asignación**

![](ElectoralSpace_files/figure-html/unnamed-chunk-20-1.png) 


**Sumas parciales en la asignación**

![](ElectoralSpace_files/figure-html/unnamed-chunk-21-1.png) 


**Historia de las electiones**

Podemos también representar datos históricos o geográficos de varias elecciones sobre el diagrama.




```r
dfvotes <- VotesData(votes=votes, election=election)

dfSpline <- generateSpline(dfvotes)

ggtern(data=df[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  #theme_rgbw()+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data = dfvotes,aes(x,y,z,color=as.factor(color)),alpha=1)+
  geom_text(data=dfvotes,aes(label=label), color="slateblue4", hjust=1.2, vjust=0.3, size=3, angle=90)+ 
  geom_text(data=df[[1]][(dots+1):(dots+nnodes),],aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  geom_path(data=dfSpline,colour="darkgreen", linetype=1, size=1)+ 
  #geom_path(data=dfvotes,colour="white", linetype=1, size=0.7)+ 
  labs(x="SocLib",y="SocCom",z="LibCon",title="Past Elections")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")
```

![](ElectoralSpace_files/figure-html/unnamed-chunk-23-1.png) 
