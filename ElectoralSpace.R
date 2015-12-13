#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)
#library(stats)
set.seed(156) #fix random generation

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
  
  index=as.matrix(1:nnodes)
  
  #node labels for seats
  label=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
  label=do.call("paste0",as.data.frame(label))
  
  
  nodes <- cbind(index,as.data.frame(nodes),label)
  names(nodes) <- c("index",letters[24:26],"label")
  return(nodes);
  
}

generateColors <- function (colorRGB, seats) {
  
  nnodes = (seats+1)*(seats+2)/2;
  nodes <- as.matrix(generateNodes(seats)[,c("x","y","z")])
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

generateDots <- function(dotsperside, method="lattice"){
  
  dots <- (dotsperside+1)*(dotsperside+2)/2
  
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
  
  else if (method == "lattice") {
    #lattice of dots
    
    #read this if you want to use dots instead of dotsperside as input variable:
    #if n is dotsperside
    #dots <- (n+1)*(n+2)/2  
    #inverse calculation:
    #n <- floor((sqrt(8*dots+1)-3)/2) #exact
    #n <- floor(sqrt(2*dots)) #aprox
    
    #then redefine number of dots
    #dots <- (n+1)*(n+2)/2
    R <- matrix(0,nrow=dots, ncol=3)
    i <- 1
    for (x in 0:dotsperside){
      
      for (y in 0:(dotsperside-x)){
        
        R[i,]<- c(x,y,dotsperside-x-y)
        i <- i+1
        
      }
    }
    R = prop.table(R,1) #rows sum 1
  }
  
  return(R);
  
}

CartesianToTernary <- function (x) {
  
  #for homogeneous simulation
  #send an subset of [0,1]x[0,1] which is an equilateral triangle to the ternary diagram 
  
  while ((x[1]<=1/2 & x[2]>2*sin(pi/3)*x[1]) | (x[1]>1/2 & x[2]>2*sin(pi/3)*(1-x[1])))
    x=runif(2);
  
  v = cbind(1-x[1]-x[2]/(2*sin(pi/3)),x[2]/sin(pi/3),x[1]-x[2]/(2*sin(pi/3)));
  
  return(v);
  
}

normalizenodes <- function (nodes) {
  
  normnodes <- nodes
  for (i in 1:dim(nodes)[1]){
    normnodes[i,] <- nodes[i,]/norm(as.matrix(nodes[i,]),"F")
  }
  
  normnodes <- as.matrix(normnodes)
  
}

UniformNearest <- function (x, seats, nodes, nnodes) {
  
  Uniform = which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,max))
  
  return(Uniform);
  
}

ManhattanNearest <- function (x, seats, nodes, nnodes) {
  
  Manhattan = which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,sum))
  
  return(Manhattan);
  
}

EuclidNearest <- function (x, seats, nodes, nnodes) {
  
  Euclid = which.min(apply((nodes/seats - rep(1,nnodes) %*% t(x))^2,1,sum))
  
  return(Euclid);
  
}

OrthodromicNearest <- function (x, seats, nodes, nnodes) {
  
  normnodes <- normalizenodes(nodes)
  Orthodromic = which.min(acos((x/norm(t(x),"F")) %*% t(as.matrix(normnodes))))
  
  return(Orthodromic);
  
}

AllocatedNode <- function (y, nodes, nnodes) {
  
  Node = which.min(apply(abs(nodes - rep(1,nnodes) %*% t(y)),1,max))
  
  return(Node);
  
}

SpatialData <- function (dotsperside, seats, step=1, threshold=0, method="lattice") {
  
  
  #Generate random electoral results
  R <- generateDots(dotsperside, method=method)
  dots=dim(R)[1]  #dots <- (dotsperside+1)*(dotsperside+2)/2
  
  Seats=max(seats)
  
  #nodes
  nodes  <- as.matrix(generateNodes(Seats)[,c("x","y","z")])
  nnodes <- (Seats+1)*(Seats+2)/2;
  
  #Indexes for Voronoi regions
  Uniform   = apply(R, 1, UniformNearest,   seats=Seats, nodes = nodes, nnodes=nnodes)
  Manhattan = apply(R, 1, ManhattanNearest, seats=Seats, nodes = nodes, nnodes=nnodes)
  Euclid    = apply(R, 1, EuclidNearest,    seats=Seats, nodes = nodes, nnodes=nnodes)
  Orthodromic=apply(R, 1, OrthodromicNearest,seats=Seats,nodes = nodes, nnodes=nnodes)
  
  
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
      
      nodes_sub <- as.matrix(generateNodes(seats[i])[,c("x","y","z")])
      
      AllocPartial[,i] = apply(S, 1, AllocatedNode, nodes=nodes_sub, nnodes=(seats[i]+1)*(seats[i]+2)/2)
      
    }
    
    #allocation order
    #computes allocation ordering for max(seats) in each step value
    AllocOrder = t(matrix(sapply(AllocStructure, function(x) as.integer(x[[j]][[length(seats)+1]])),nrow=Seats,ncol=dots))
    vectorOrderCode=as.matrix(3^(0:(Seats-1)))
    #this coding highlights the last seats over the first ones
    # so that there is contrast between adjacent regions
    AllocOrderCode = (matrix(AllocOrder, ncol=Seats)-1) %*% as.matrix(vectorOrderCode)
    
    dfPartial=as.data.frame(AllocPartial)
    names(dfPartial) <- do.call(paste0,as.data.frame(cbind("All",seats)))
    names(dfPartial)[which(seats==Seats)] <- "Allocated"
    
    #assembling a data frame
    df0 = data.frame(
      
      
      x             = as.matrix(R[,1]),
      y             = as.matrix(R[,2]),
      z             = as.matrix(R[,3]),
      
      Sx            = S[,1],
      Sy            = S[,2],
      Sz            = S[,3],
      
      Euclid,
      Manhattan,
      Uniform,
      Orthodromic,
      
      Malapportionment = AllocPartial[1:dots,length(seats)] != Euclid[1:dots],
      Malapportionment2 = AllocPartial[1:dots,length(seats)] != Orthodromic[1:dots],
      
      AllocOrderCode
      
      
    )
    
    df0=cbind(df0,dfPartial)
    
    df[[j]]<-df0
    
  } #end of loop for different 'step's
  
  
  return(df);
  
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

VotesData <- function (votes=matrix(0,1), election=matrix(0,1)) {
  
  votes <- prop.table(votes,1)  
  
  x <- votes[,1]
  y <- votes[,2]
  z <- votes[,3]
  
  el=dim(votes)[1]
  label = as.character(election)
  
  
  #assembling a data frame
  df = data.frame(
    
    type = matrix("vote",el,1),
    
    x ,
    y ,
    z ,
    
    label,
    
    election,
    
    color = matrix(NA,el,1)
    
  )
  
  return(df);
  
}

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

modifiedlog <- function(x) {
  
  #modified log function for the calculation of entropy
  #avoids log(0)=NaN
  
  x <- as.matrix(x)
  r<-log(x)
  r[is.finite(r)==F] <-0
  return(r)
  
}


#plot example
ggtern(data=data.frame(x=0.2,y=0.3,z=0.5),aes(x,y,z))+
  theme_arrowdefault()+geom_point()

#colors for palettes
colorRGB1 <- c(242,74,87, 124,218,198, 91,168,246) #c("#f24a57","#7cdac6","#5ba8f6")
colorRGB2 <- c(226,10,23, 100,207,151, 0,112,184) #c("#e20a17","#64cf97","#0070b8")
colorRGB3 <- c(242,74,87, 124,218,198, 2,114,195) #<--
colorRGB4 <- c(256,0,0, 0,256,0, 0,0,256)
colorRGB5 <- c(224,0,0, 96,224,176, 96,128,224)
colorRGB5 <- c(240,0,0, 96,224,176, 96,128,224)
colorRGB6 <- c(256,256,0, 0,256,256, 256,0,256) ##YMC test shifted colours
colorRGB7 <- c(256,256,256, 256,0,0, 0,0,0) #WRB

colorRGB0 <- colorRGB3 #choose color palette

seats=5
NodesData <- generateNodes(max(seats))

#plot nodes for seat allocation
ggtern(data=NodesData,aes(x,y,z,color=as.factor(index)))+
  theme_minimal()+
  geom_point(alpha=1, size=5)+
  geom_text(aes(label=label,color=as.factor(index)), hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Nodes for seat allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")


#Allocation example (step=2 Sainte-Lagu?; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
#alloc(letters[1:3], votes, seats=5, step=1)
#alloc(letters[1:3], votes, seats=c(3,5,4), step=c(2,1), threshold=0.5)
alloc(letters[1:3], votes, 9, c(1,2), .05) #print seats sum and allocation


#presets
seats=2:5;
step=c(1,2); #(2 Sainte-Lagu? 1 D'Hondt)
dotsperside=199 #dots <- (dotsperside+1)*(dotsperside+2)/2
threshold=0

#Spatial data

dots <- (dotsperside+1)*(dotsperside+2)/2

df = SpatialData(dotsperside, seats, step)

#df = SpatialData(dotsperside, 5, 1, threshold)

#df = SpatialData(dotsperside, c(3,5,4), c(2,1), threshold)

dfT = SpatialData(dotsperside, seats=5, threshold=.20)


head(df[[1]][sample(1:dots,10,replace=F),]) #sample data for step=1 and seats=5



#Allocation

a1 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="D'Hondt Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")

a2 <- ggtern(data=df[[2]],aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")

ggtern.multi(a1, a2, cols=2)

#Voronoi

v1 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(Manhattan)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Voronoi Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")

#Orthodromic Voronoi

v2 <- ggtern(data=df[[2]],aes(x,y,z,color=as.factor(Orthodromic)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Orthodromic Voronoi Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")

ggtern.multi(a2, v1, cols=2)
ggtern.multi(a2, v2, cols=2)


#equivalence between Voronoi regions using different distances and between those and electoral regions
sum(df[[1]]$Euclid[1:dots]==df[[1]]$Manhattan[1:dots])/dots
sum(df[[1]]$Uniform[1:dots]==df[[1]]$Manhattan[1:dots])/dots
sum(df[[1]]$Uniform[1:dots]==df[[1]]$Euclid[1:dots])/dots

#Differences between planar and orthodromic Voronoi regions
sum(df[[1]]$Euclid[1:dots]==df[[1]]$Orthodromic[1:dots])/dots

ggtern(data=df[[1]],aes(x,y,z,color=Malapportionment!=Malapportionment2)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey30")+
  labs(x="X",y="Y",z="Z",title="Planar vs Orthodromic regions")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "Greys", na.value = "grey60", guide = FALSE)


#size of the regions # 1 ~= dots/nnodes
nnodes=(max(seats)+1)*(max(seats)+2)/2;
RegionSize  = table(df[[1]]$Allocated[1:dots])/(dots/nnodes)
RegionSize2 = table(df[[2]]$Allocated[1:dots])/(dots/nnodes)

par(mfrow=c(2,1))
plot(RegionSize, main="D'Hondt region sizes")
plot(RegionSize2,main="Sainte-Laguë region sizes" )
par(mfrow=c(1,1))


#malapportionment
m1 <- ggtern(data=df[[1]],aes(x,y,z,color=Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="D'Hondt Malapportionment")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

m2 <- ggtern(data=df[[2]],aes(x,y,z,color=Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Malapportionment")+
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

ggtern.multi(m1, m2, cols=2)

#orthodromic malapportionment
om1 <- ggtern(data=df[[1]],aes(x,y,z,color=Malapportionment2)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="D'Hondt Malapportionment")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

om2 <- ggtern(data=df[[2]],aes(x,y,z,color=Malapportionment2)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Malapportionment")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

ggtern.multi(om1, om2, cols=2)

#points not allocated by D'Hondt in their corresponding Voronoi region
sum(df[[1]]$Malapportionment[1:dots])/dots
#points not allocated by Sainte-Lagu? in their corresponding Voronoi region
sum(df[[2]]$Malapportionment[1:dots])/dots

#points not allocated by D'Hondt in their corresponding orthodromic Voronoi region
sum(df[[1]]$Malapportionment2[1:dots])/dots
#points not allocated by Sainte-Laguë in their corresponding orthodromic Voronoi region
sum(df[[2]]$Malapportionment2[1:dots])/dots


#Threshold effect
t1 <- ggtern(data=dfT[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Threshold effect on regions")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")

ggtern.multi(a1, t1, cols=2)


#ordering subregions
o1 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(AllocOrderCode))) +
  theme_bw()+
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),color="black")+
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  scale_colour_grey(start = 0.1, end = 1, na.value = "black", guide = FALSE)

o2 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(AllocOrderCode))) +
  theme_bw()+
  geom_point(alpha=1) +
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  labs(x="X",y="Y",z="Z",title="They come in colors")+
  scale_colour_manual(values=
      generateOrderColors(colorRGB0,max(seats),sort(unique(df[[1]]$AllocOrderCode))),
      guide=FALSE, na.value="khaki2")

ggtern.multi(o1, o2, cols=2)

#Partial Allocations

p1 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(df[[1]]$All2)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 2 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,2), guide=FALSE)

p2 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(df[[1]]$All3)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 3 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,3), guide=FALSE)

p3 <- ggtern(data=df[[1]],aes(x,y,z,color=as.factor(df[[1]]$All4)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 4 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,4), guide=FALSE)

p4 <- ggtern(data=df[[2]],aes(x,y,z,color=as.factor(df[[2]]$All2)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Lagu?, 2 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,2), guide=FALSE)

p5 <- ggtern(data=df[[2]],aes(x,y,z,color=as.factor(df[[2]]$All3)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Lagu?, 3 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,3), guide=FALSE)

p6 <- ggtern(data=df[[2]],aes(x,y,z,color=as.factor(df[[2]]$All4)))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Lagu?, 4 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,4), guide=FALSE)

plist=list(p1,p4,p2,p5,p3,p6)
ggtern.multi(plotlist=plist, cols=3)



#history of election results


#nelect=10;
#votes=matrix(runif(nelect*3),nelect,3)
#votes <- prop.table(votes,1)

election=as.matrix(seq(from=1979, to=2015, by=4))

#example
#seats=5; step=1;
nelect=10;
votes=matrix(0,nelect,3)
votes[,1]=c(350,430,290,400,390,450,470,360,280,230)
votes[,3]=c(250,180,220,160,140,110, 80, 60, 90,140)
votes[,2]=c( 80, 90,130, 70, 90, 80,120,140,130,110)

#votes[,1]=c(355,437,282,393,408,479,468,356,287,234)
#votes[,3]=c(258,  0,185,144,159, 90, 47, 73,114, 75)
#votes[,2]=c(  0,  0, 82,  0,  0,  0,  0,  0,  0,138)

dfvotes <- VotesData(votes=votes, election=election)

dfSpline <- generateSpline(dfvotes)

ggtern(data=df[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  #theme_rgbw()+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_point(data=dfvotes,aes(x,y,z),color="orange",alpha=1)+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  geom_text(data=dfvotes,aes(label=label), color="slateblue4", hjust=1.2, vjust=0.3, size=3, angle=90)+ 
  geom_path(data=dfSpline,colour="darkgreen", linetype=1, size=1)+ 
  #geom_path(data=dfvotes,colour="white", linetype=1, size=0.7)+ 
  labs(x="SocLib",y="SocCom",z="LibCon",title="Past Elections")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")



#other colour palettes
library(RColorBrewer)
# RdYlGn 11
# Spectral 11
# YlGnBu 9
values=colorRampPalette(brewer.pal(8,"Spectral"))(nnodes)
values=rev(colorRampPalette(brewer.pal(8,"PuBuGn"))(nnodes))
values=colorRampPalette(brewer.pal(9,"Spectral"))(nnodes+10)[14:(nnodes+13)]
values = colorRampPalette(as.character(c("#00FF00","#FF0000","#0000FF")))(nnodes)
values = colorRampPalette(c(rev(brewer.pal(9, "Blues")),brewer.pal(9, "Reds")))(nnodes)
values = colorRampPalette(as.character(c("#f24a57","#7cdac6","#5ba8f6")))(nnodes)
values = colorRampPalette(as.character(c("#e20a17","#64cf97","#0070b8")))(nnodes)


#install.packages("nnet")
library('nnet')
fit <- multinom(Allocated ~ x + y + z, data = df)

summary(fit)
str(fit)
plot(fit, df)
print(fit)
fit$coefs
fit$terms
fit$fitted


#Entropy
Seats <- length(seats)
df[[1]][,11+Seats+1] <- rowSums(-df[[1]][,1:3]*modifiedlog(df[[1]][,1:3]))
df[[1]][,11+Seats+2] <- rowSums(-df[[1]][,4:6]/max(seats)*modifiedlog(df[[1]][,4:6]/max(seats)))

#Effective number of parties
df[[1]][,11+Seats+3] <- 1/rowSums((df[[1]][,1:3])^2)
df[[1]][,11+Seats+4] <- 1/rowSums((df[[1]][,4:6]/max(seats))^2)

names(df[[1]])[(11+Seats+1):(11+Seats+4)] <- c('Entropy','CameraEntropy','Parties','CameraParties')

options(warn=-1) 

e1 <- ggtern(data=df[[1]],aes(x,y,z,color=Entropy))+
  theme_bw()+
  geom_point(alpha=1)+
  #geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  #geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Votes Entropy")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  #scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")
  scale_colour_gradient2(low='brown3', mid="aquamarine3", high='black', midpoint=log(2), na.value = "grey50", guide = "colourbar")

e2 <- ggtern(data=df[[1]],aes(x,y,z,color=CameraEntropy))+
  theme_bw()+
  geom_point(alpha=1)+
  #geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  #geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Camera Entropy")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  #scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")
  scale_colour_gradient2(low='brown3', mid="aquamarine3", high='black', midpoint=log(2), na.value = "grey50", guide = "colourbar")

ggtern.multi(e1, e2, cols=2)


n1 <- ggtern(data=df[[1]],aes(x,y,z,color=Parties))+
  theme_bw()+
  geom_point(alpha=1)+
  #geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  #geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Eff. num of parties")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  #scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")
  scale_colour_gradient2(low='brown3', mid="aquamarine3", high='black', midpoint=2, na.value = "grey50", guide = "colourbar")

n2 <- ggtern(data=df[[1]],aes(x,y,z,color=CameraParties))+
  theme_bw()+
  geom_point(alpha=1)+
  #geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  #geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Eff. num of parties")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  #scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")
  scale_colour_gradient2(low='brown3', mid="aquamarine3", high='black', midpoint=2, na.value = "grey50", guide = "colourbar")

ggtern.multi(n1, n2, cols=2)

ggtern.multi(e1, e2, n1, n2, cols=2)



#Largest remainder method with Hare Quota is Voronoi

alloc <- function(parties, votes, seats, step=1, threshold=0){
  
  #function for Hare Quota seat allocation
  votes=votes*(votes>=(threshold*sum(votes)))
  
  nparties = length(votes)
  
  quota = sum(votes)/seats
  QAllocation = floor(votes/quota)
  rests = votes %% quota
  
  qtable <- data.frame(
    parties = 1:nparties,
    rests,
    votes
  )
  
  if (seats>sum(QAllocation)){
    
    mayorRestsParties <- qtable[order(-rests,-votes),]$parties[1:(seats - sum(QAllocation))]
    mayorRestsLogic <- is.element(1:nparties, mayorRestsParties)
    
  } else {mayorRestsLogic <- matrix(F,1,nparties)}
  
  
  Allocation = list(list(QAllocation + mayorRestsLogic, matrix(parties[1],seats)))
  #ordering is not taken into account
  
  return(Allocation);
  
}

seats=5;

dfH = SpatialData(dotsperside, seats)

#head(dfH[[1]][sample(1:dots,10,replace=F),])

#NodesData <- generateNodes(max(seats))

ggtern(data=dfH[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Largest remainder - Hare Quota Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,max(seats)), guide=FALSE, na.value="khaki2")


#malapportionment
mH1 <- ggtern(data=dfH[[1]],aes(x,y,z,color=Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="Hare Malapportionment")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

#orthodromic malapportionment
mH2 <- ggtern(data=dfH[[1]],aes(x,y,z,color=Malapportionment2)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_point(data=NodesData,aes(x,y,z),alpha=0.8,color="grey60")+
  labs(x="X",y="Y",z="Z",title="Hare Malapportionment (Orthodromic)")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

ggtern.multi(mH1, mH2, cols=2)

#points not allocated by Hare in their corresponding Voronoi region
sum(dfH[[1]]$Malapportionment[1:dots])/dots
#points not allocated by Hare in their corresponding orthodromic Voronoi region
sum(dfH[[1]]$Malapportionment2[1:dots])/dots






