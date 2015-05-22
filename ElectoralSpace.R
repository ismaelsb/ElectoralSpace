#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)
#library(stats)

alloc <- function(parties, votes, seats, step, threshold, partial=seats){
  
  #function for the seat allocation and its ordering
  votes=votes*(votes>=(threshold*sum(votes)))
  quotienstable <- data.frame( 
    parties    = rep(parties, each = seats), 
    quotients  = as.vector(sapply(votes, function(x) 
      x/seq(from=1, to=1+step*(seats-1), by=step) )),
    votesrep   = rep(votes, each = seats)
  ) 
  quotienstable <- quotienstable$parties[order(-quotienstable$quotients, -quotienstable$votesrep)] [1:seats] 
  
  SeatsList=list()
  
  #partial=TRUE (default partial=seats if don't want to compute the partial sums)
  for (i in partial:seats){
    SeatsList[[i]]= table(quotienstable[1:i]) #all seats allocated for i=seats 
  }
  
  SeatsList[[seats+1]] = as.matrix(quotienstable) #ordering
  
  SeatsList = SeatsList[partial: (seats+1)]
  
  return(SeatsList);
}

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
    R = R/rowSums(R)
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

UniformNearest <- function (x, nodes, nnodes) {
  
  Uniform = which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,max))
  
  return(Uniform);
  
}

ManhattanNearest <- function (x, nodes, nnodes) {
  
  Manhattan = which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,sum))
  
  return(Manhattan);
  
}

EuclidNearest <- function (x, nodes, nnodes) {
  
  Euclid = which.min(apply((nodes/seats - rep(1,nnodes) %*% t(x))^2,1,sum))
  
  return(Euclid);
  
}

AllocatedNode <- function (y, nodes, nnodes) {
  
  Node = which.min(apply(abs(nodes - rep(1,nnodes) %*% t(y)),1,max))
  
  return(Node);
  
}

generateSpline <- function (dfvotes, dots, seats, nelect, method="natural") {
  
  #generates a spline curve through elections history
  nnodes = (seats+1)*(seats+2)/2;
  
  splineX <- spline(x=1:nelect,y=dfvotes[(dots+nnodes+1):(dots+nnodes+nelect),]$x, method=method)
  splineY <- spline(x=1:nelect,y=dfvotes[(dots+nnodes+1):(dots+nnodes+nelect),]$y, method=method)
  
  dfSpline = matrix(0, length(splineX[[1]]), 3)
  dfSpline[,1] <- splineX[[2]]
  dfSpline[,2] <- splineY[[2]]
  dfSpline[,3] <- 1-dfSpline[,1]-dfSpline[,2]
  dfSpline = as.data.frame(dfSpline)
  names(dfSpline) <- letters[24:26]
  
  return(dfSpline);
  
}

AssembleData <- function (dots, seats, step, threshold, vectorAllocCode = c((seats+1)^2,seats+1,1), vectorOrderCode=as.matrix(3^(0:(seats-1))), votes=matrix(0,1), election=matrix(0,1), partial=seats) {
  
  #Generate random electoral results
  R <- generateDots(dots)
  
  #nodes
  nodes  = generateNodes(seats)
  nnodes = (seats+1)*(seats+2)/2;
  
  #allocate seats
  AllocStructure <- apply(R, 1, function(x) alloc(as.character(1:3), x, seats, step, threshold, partial=partial))
  #input partial=1 to compute partial sums or partial=seats to compute only seats
  
  #allocation
  AllocPartial=matrix(0,dots,seats-partial+1)
  
  for (i in partial:seats){
    S = t(matrix(sapply(AllocStructure, function(x) x[[i-partial+1]]),nrow=3,ncol=dots))
    vectorAllocCode = c((i+1)^2,i+1,1)
    AllocPartial[,i-partial+1] = matrix(S, ncol=3) %*% as.matrix(vectorAllocCode)
  }
  
  AllocCode=AllocPartial[,seats-partial+1]
  
  #allocation order
  AllocOrder = t(matrix(sapply(AllocStructure, function(x) as.integer(x[[seats-partial+2]])),nrow=seats,ncol=dots))
  AllocOrderCode = (matrix(AllocOrder, ncol=seats)-1) %*% as.matrix(vectorOrderCode)
  
  #Indexes for Voronoi and Allocation regions
  Uniform   = apply(R, 1, UniformNearest,   nodes = nodes, nnodes=nnodes)
  Manhattan = apply(R, 1, ManhattanNearest, nodes = nodes, nnodes=nnodes)
  Euclid    = apply(R, 1, EuclidNearest,    nodes = nodes, nnodes=nnodes)
  Allocated = apply(S, 1, AllocatedNode,    nodes = nodes, nnodes=nnodes)
  
  #node labels for seats
  label=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
  label=do.call("paste0",as.data.frame(label))
  label=rbind(matrix("",dots),as.matrix(label))
  
  #elections
  elect = matrix(NA,dots+nnodes)
  if (norm(election)!=0) elect = rbind(elect, election)
  
  #attach node points                
  x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/seats))
  y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/seats))
  z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/seats))
  
  el=0;
  if (norm(votes) != 0){
    
    x <- c(x, votes[,1])
    y <- c(y, votes[,2])
    z <- c(z, votes[,3])
    el=dim(votes)[1]
    label = c(label,as.character(election))
    
  }
  
  AllocPartial=rbind(AllocPartial,matrix(NA,nnodes+el,seats-partial+1))
  dfPartial=as.data.frame(AllocPartial)
  
  
  #assembling a data frame
  df = data.frame(
    
    type          = c(matrix("dot",dots), matrix("node",nnodes), matrix("vote",el,1)),
    
    x ,
    y ,
    z ,
    
    AllocCode     =  c(AllocCode, matrix(NA,nnodes+el)),
    #code          = c(rgb(S/seats), matrix(NA,nnodes+el)),
    
    AllocOrderCode=  c(AllocOrderCode, matrix(NA,nnodes+el)),
    #AllocOrderCodeC= c(rgb(AllocOrderCode/(3^seats),0,0), matrix(NA,nnodes+el)),
    
    Euclid        = c(Euclid, matrix(NA,nnodes+el)),
    #codeEuclid    = c(rgb(nodes[Euclid[1:dots],]/seats), matrix(NA,nnodes+el)),
    
    Manhattan     = c(Manhattan, matrix(NA,nnodes+el)),
    #codeManhattan = c(rgb(nodes[Manhattan[1:dots],]/seats), matrix(NA,nnodes+el)),
    
    Uniform       = c(Uniform, matrix(NA,nnodes+el)),
    #codeUniform   = c(rgb(nodes[Uniform[1:dots],]/seats), matrix(NA,nnodes+el)),
    
    Allocated     = c(Allocated, matrix(NA,nnodes+el)),
    Malapportionment = c((Allocated!=Euclid)[1:dots], matrix(NA,nnodes+el)),
    
    label,
    
    elect
  )
  
  df=cbind(df,dfPartial)
  
  return(df);
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
alloc(letters[1:3], votes, 10, 1, .05) #print seats sum and allocation

#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=50000
threshold=.05

#Assemble the data frame
#df = AssembleData(dots, seats, step, threshold)
df = AssembleData(dots, seats, step, threshold, partial=TRUE)

df2 =AssembleData(dots, seats, step=2,threshold, partial=2)

#alternative encodings for Alloc and AllocOrder:
#default encoding:
#vectorAllocCode = c((seats+1)^2,seats+1,1)
#vectorOrderCoder = as.matrix(3^(0:(seats-1)))
#reverse coding for allocation
#vectorAllocCoderev = c(1,seats+1,(seats+1)^2)
#reverse encoding for ordering regions
#vectorOrderCoderev = as.matrix(rev(3^(0:(seats-1))))
#df = AssembleData(dots, seats, step, threshold, vectorAllocCode=vectorAllocCoderev)
#df = AssembleData(dots, seats, step, threshold, vectorOrderCode=vectorOrderCoderev)

#Allocation
ggtern(data=df,aes(x,y,z,color=as.factor(AllocCode))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="D'Hondt Allocation")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

ggtern(data=df2,aes(x,y,z,color=as.factor(AllocCode))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Allocation")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

#Voronoi
ggtern(data=df,aes(x,y,z,color=as.factor(Manhattan))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Voronoi Allocation")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

#equivalence between Voronoi regions using different distances and between those and electoral regions
sum(df$Euclid[1:dots]==df$Manhattan[1:dots])/dots
sum(df$Uniform[1:dots]==df$Manhattan[1:dots])/dots
sum(df$Uniform[1:dots]==df$Euclid[1:dots])/dots

#points not allocated by D'Hondt in their corresponding Voronoi region
sum(df$Malapportionment[1:dots])/dots
#points not allocated by Sainte-Laguë in their corresponding Voronoi region
sum(df2$Malapportionment[1:dots])/dots

#size of the regions # 1 ~= dots/nnodes
nnodes=(seats+1)*(seats+2)/2;
RegionSize = table(df$Allocated[1:dots])/(dots/nnodes)
RegionSize2 = table(df2$Allocated[1:dots])/(dots/nnodes)

par(mfrow=c(2,1))
plot(RegionSize, main="D'Hondt region sizes")
plot(RegionSize2,main="Sainte-Laguë region sizes" )
par(mfrow=c(1,1))

#malapportionment
m1 <- ggtern(data=df,aes(x,y,z,color=!Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="D'Hondt Malapportionment")+
  scale_colour_grey(na.value = "black", guide = FALSE)

m2 <- ggtern(data=df2,aes(x,y,z,color=!Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Malapportionment")+
  scale_colour_grey(na.value = "black", guide = FALSE)

ggtern.multi(m1, m2, cols=2)

#ordering subregions
ggtern(data=df,aes(x,y,z,color=as.factor(AllocOrderCode))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  scale_colour_grey(start = 0.1, end = 1, na.value = "black", guide = FALSE)

#Partial Allocations

p1 <-  ggtern(data=df[1:dots,],aes(x,y,z,color=as.factor(df[1:dots,13+2]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 2 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

p2 <-  ggtern(data=df[1:dots,],aes(x,y,z,color=as.factor(df[1:dots,13+3]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 3 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

p3 <-  ggtern(data=df[1:dots,],aes(x,y,z,color=as.factor(df[1:dots,13+4]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 4 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

p4 <-  ggtern(data=df2[1:dots,],aes(x,y,z,color=as.factor(df2[1:dots,13+1]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 2 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

p5 <-  ggtern(data=df2[1:dots,],aes(x,y,z,color=as.factor(df2[1:dots,13+2]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 3 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

p6 <-  ggtern(data=df2[1:dots,],aes(x,y,z,color=as.factor(df2[1:dots,13+3]))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 4 seats")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

plist=list(p1,p4,p2,p5,p3,p6)
ggtern.multi(plotlist=plist, cols=3)

#history of election results
#example
#seats=7; step=1;
nelect=10;
votes=matrix(runif(nelect*3),nelect,3)
votes <- votes/rowSums(votes)
election=as.matrix(seq(from=1979, to=2015, by=4))

#example
seats=7; step=1;
nelect=10;
votes=matrix(0,nelect,3)
votes[,1]=c(350,430,290,400,390,450,470,360,280,230)
votes[,3]=c(250,180,220,160,140,110, 80, 60, 90,140)
votes[,2]=c( 80, 90,130, 70, 90, 80,120,140,130,110)
votes <- votes/rowSums(votes)

dfvotes = AssembleData(dots, seats, step, threshold, votes=votes, election=election)

dfSpline <- generateSpline(dfvotes, dots, seats, nelect)

ggtern(data=dfvotes,aes(x,y,z,color=as.factor(AllocCode))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  geom_text(aes(label=label), color="red", hjust=0.5,vjust=-0.6, size=3)+ 
  geom_path(data=dfSpline,colour="darkgreen", linetype=1, size=1)+  
  labs(x="SocLib",y="SocCom",z="LibCon",title="Past Elections")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)

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


