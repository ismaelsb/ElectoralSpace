#Uncommmet to install ggtern from an online CRAN repository:
#install.packages("ggtern")
#install.packages("rgl")
#Load the ggtern library:
library(ggtern)
#library(stats)
library("rgl")
set.seed(156) #fix random generation


#Transform data from 4-simplex to 3D space
SimplexCoordinates <- c(0, 0, 0,   1, 0, 0,   1/2, sqrt(3)/2, 0,   1/2, sqrt(3)/6, sqrt(2/3))
SimplexToCartesian <- matrix(SimplexCoordinates,nrow=4,ncol=3,byrow=T)



generate4Nodes <- function(seats){ 
  
  #define nodes
  nnodes = (seats+1)*(seats+2)*(seats+3)/6;
  nodes  = matrix(0,nnodes,4);
  t=1;
  for (i in 0:seats){
    for (j in i:seats){
      for (k in j:seats){
        nodes[t,] = c(i,j-i,k-j,seats-k);
        t=t+1;
      }
    }
  }
  
  normnodes = prop.table(nodes,1) #rows sum 1
  NodesCoords <- normnodes %*% SimplexToCartesian
  
  index=as.matrix(1:nnodes)
  
  #node labels for seats
  label=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]),matrix("-",nnodes),as.character(nodes[,4]))
  label=do.call("paste0",as.data.frame(label))
  
  
  nodes <- cbind(index, NodesCoords, as.data.frame(nodes),label)
  names(nodes) <- c("index", letters[24:26], letters[1:4],"label")
  return(nodes);
  
}

generateColors4 <- function (colorRGB, seats) {
  
  nnodes = (seats+1)*(seats+2)*(seats+3)/6;
  nodes <- as.matrix(generateNodes(seats)[,c("a","b","c","d")])
  #decimal codes for colors in nodes
  Code4 = floor(nodes*255.9/seats) #255.9 avoids seats -> 256 -> HEX #100 case
  #linear transformation to fixed extreme colors
  CodeDecRGB = floor(Code4%*%matrix(colorRGB,4,3, byrow=T)/256)
  #hex codes for colors in nodes
  CodeRGB=cbind(matrix("#",nnodes),format(as.hexmode(CodeDecRGB[,1]),width=2),format(as.hexmode(CodeDecRGB[,2]),width=2),format(as.hexmode(CodeDecRGB[,3]),width=2))
  CodeRGB=do.call("paste0",as.data.frame(CodeRGB))
  values=CodeRGB
  
  return(values);
  
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
  
  else if (method == "tetrahedron") {
    
    dots <- (dotsperside+1)*(dotsperside+2)*(dotsperside+3)/6
    
    R <- matrix(0,nrow=dots, ncol=4)
    i <- 1
    for (w in 0:dotsperside){
      
      for (x in 0:(dotsperside-w)){
        
        for (y in 0:(dotsperside-w-x)){
          
          R[i,]<- c(w,x,y,dotsperside-w-x-y)
          i <- i+1
          
          
        }
      }
    }
    R = prop.table(R,1) #rows sum 1
    
    
  }
  
  return(R);
  
}

SpatialData4 <- function (dotsperside, seats, step=1, threshold=0, method="tetrahedron") {
  
  
  #Generate random electoral results
  R <- generateDots(dotsperside, method=method)
  dots=dim(R)[1]  #dots <- (dotsperside+1)*(dotsperside+2)*(dotsperside+3)/6
  
  Coords <- R %*% SimplexToCartesian
  
  
  Seats=max(seats)
  
  #nodes
  nodes  <- as.matrix(generateNodes(Seats)[,c("a","b","c","d")])
  nnodes <- (Seats+1)*(Seats+2)*(Seats+3)/6;
  
  #Indexes for Voronoi regions
  Uniform   = apply(R, 1, UniformNearest,   seats=Seats, nodes = nodes, nnodes=nnodes)
  Manhattan = apply(R, 1, ManhattanNearest, seats=Seats, nodes = nodes, nnodes=nnodes)
  Euclid    = apply(R, 1, EuclidNearest,    seats=Seats, nodes = nodes, nnodes=nnodes)
  Orthodromic=apply(R, 1, OrthodromicNearest,seats=Seats,nodes = nodes, nnodes=nnodes)
  
  
  #allocate seats
  AllocStructure <- apply(R, 1, function(x) alloc(as.character(1:4), x, seats, step, threshold))
  #input a vector of values for seats of to compute partial sums
  
  df=list() #list of dataframes, one for each 'step' value
  
  #loop for diferent 'step's
  for (j in 1:length(step)){
    
    #allocation
    AllocPartial=matrix(0,dots,length(seats))
    
    for (i in 1:length(seats)){
      
      S = t(matrix(sapply(AllocStructure, function(x) x[[j]][[i]]),nrow=4,ncol=dots))
      
      nodes_sub <- as.matrix(generateNodes(seats[i])[,c("a","b","c","d")])
      
      AllocPartial[,i] = apply(S, 1, AllocatedNode, nodes=nodes_sub, nnodes=(seats[i]+1)*(seats[i]+2)*(seats[i]+3)/6)
      
    }
    
    #allocation order
    #computes allocation ordering for max(seats) in each step value
    AllocOrder = t(matrix(sapply(AllocStructure, function(x) as.integer(x[[j]][[length(seats)+1]])),nrow=Seats,ncol=dots))
    vectorOrderCode=as.matrix(4^(0:(Seats-1)))
    #this coding highlights the last seats over the first ones
    # so that there is contrast between adjacent regions
    AllocOrderCode = (matrix(AllocOrder, ncol=Seats)-1) %*% as.matrix(vectorOrderCode)
    
    dfPartial=as.data.frame(AllocPartial)
    names(dfPartial) <- do.call(paste0,as.data.frame(cbind("All",seats)))
    names(dfPartial)[which(seats==Seats)] <- "Allocated"
    
    #assembling a data frame
    df0 = data.frame(
      
      x             = as.matrix(Coords[,1]),
      y             = as.matrix(Coords[,2]),
      z             = as.matrix(Coords[,3]),
      
      x1             = as.matrix(R[,1]),
      x2             = as.matrix(R[,2]),
      x3             = as.matrix(R[,3]),
      x4             = as.matrix(R[,4]),
      
      S1            = S[,1],
      S2            = S[,2],
      S3            = S[,3],
      S4            = S[,3],
      
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

generate4Colors <- function (color4RGB, seats, df) {
  
  Seats <- max(seats)
  ColorMatrix = matrix(color4RGB, 4,3, byrow=T)
  NodesData <- generate4Nodes(Seats)
  
  dColors <- as.data.frame(as.matrix(NodesData[df$Allocated,][,5:8]/Seats)%*%ColorMatrix/256)
  names(dColors) <- c("R","G","B")
  
  return(dColors);
  
}

#Voronoi
generate4VColors <- function (color4RGB, seats, df) {
  
  Seats <- max(seats)
  ColorMatrix = matrix(color4RGB, 4,3, byrow=T)
  NodesData <- generate4Nodes(Seats)
  
  dColors <- as.data.frame(as.matrix(NodesData[df$Uniform,][,5:8]/Seats)%*%ColorMatrix/256)
  names(dColors) <- c("R","G","B")
  
  return(dColors);
  
}



seats=5
NodesData <- generate4Nodes(max(seats))

#colors for palettes
color4RGB <- c(242,74,87, 122,55,139,  2,114,195, 255,127,0) #<--RPBO
color4RGB <- c(242,74,87, 124,218,198, 2,114,195, 130,20,070) #<--RGBP
color4RGB <- c(242,74,87, 124,218,198, 2,114,195, 210,230,10) #<--RGBY

#plot nodes for seat allocation
ggplot(data=NodesData,aes(x,y,z,color=as.factor(index)))+
  theme_minimal()+
  geom_point(alpha=1, size=5)+
  geom_text(aes(label=label,color=as.factor(index)), hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Nodes for seat allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors4(color4RGB,max(seats)), guide=FALSE, na.value="khaki2")


#presets
seats=2:5;
step=c(1,2); #(2 Sainte-Laguë 1 D'Hondt)
dotsperside=64 #dots <- (dotsperside+1)*(dotsperside+2)*(dotsperside+3)/6
threshold=0

df4 = SpatialData4(dotsperside, seats, step)



ggplot(data=df4[[1]],aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=.5)+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="D'Hondt Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors4(color4RGB,max(seats)), guide=FALSE, na.value="khaki2")


#size of the regions # 1 ~= dots/nnodes
nnodes=(max(seats)+1)*(max(seats)+2)*(max(seats)+2)/6
dots <- (dotsperside+1)*(dotsperside+2)*(dotsperside+3)/6
RegionSize  = table(df4[[1]]$Allocated[1:dots])/(dots/nnodes)
RegionSize2 = table(df4[[2]]$Allocated[1:dots])/(dots/nnodes)

par(mfrow=c(2,1))
plot(RegionSize, main="D'Hondt region sizes")
plot(RegionSize2,main="Sainte-Laguë region sizes" )
par(mfrow=c(1,1))



#Plots for D'Hondt, Sainte-Laguë and Voronoi regions
open3d()

plot3d(df4[[1]]$x, df4[[1]]$y, df4[[1]]$z, 
       col=rgb(generate4Colors(color4RGB, seats, df4[[1]])),
       xlab="", ylab="", zlab="",
       size=3, box=F, axes=F, top=T)
play3d( spin3d(axis=c(2,2,3), rpm=15), duration = 5 )


plot3d(df4[[1]]$x, df4[[1]]$y, df4[[1]]$z, 
       col=rgb(generate4Colors(color4RGB, seats, df4[[2]])),
       xlab="", ylab="", zlab="",
       size=3, box=F, axes=F, top=T)
play3d( spin3d(axis=c(2,2,3), rpm=15), duration = 5 )


plot3d(df4[[1]]$x, df4[[1]]$y, df4[[1]]$z, 
       col=rgb(generate4VColors(color4RGB, seats, df4[[1]])),
       xlab="", ylab="", zlab="",
       size=3, box=F, axes=F, top=T)
play3d( spin3d(axis=c(2,2,3), rpm=15), duration = 5 )

