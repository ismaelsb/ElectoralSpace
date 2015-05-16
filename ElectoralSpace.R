#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)

#defining a function for the seat allocation
alloc <- function(parties, votes, seats, step, threshold){ 
  votes=votes*(votes>=(threshold*sum(votes)))
  qtable <- data.frame( 
    parties = rep(parties, each = seats), 
    quotients  = as.vector(sapply(votes, function(x) x / 
                                    seq(from=1, to=1+step*(seats-1), by=step) )),
    votesrep = rep(votes, each = seats)
  ) 
  qtable <- qtable$parties[order(-qtable$quotients, -qtable$votesrep)] [1:seats] 
  table(qtable) 
}

#defining a function for the seat allocation order
allocOrder <- function(parties, votes, seats, step, threshold){ 
  votes=votes*(votes>=(threshold*sum(votes)))
  qtable <- data.frame( 
    parties = rep(parties, each = seats), 
    quotients  = as.vector(sapply(votes, function(x) x / 
                                    seq(from=1, to=1+step*(seats-1), by=step) )),
    votesrep = rep(votes, each = seats)
  )
  
  qtable <- qtable$parties[order(-qtable$quotients, -qtable$votesrep)] [1:seats] 
  #table(qtable)
  as.matrix(qtable)  
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

AssembleData <- function (R, seats, step, threshold ,votes=matrix(0,1), dates=matrix(0,1)) {
  
  dots=dim(R)[1]
  #define nodes
  nnodes=(seats+1)*(seats+2)/2;
  nodes= matrix(0,nnodes,3);
  t=1;
  for (i in 0:seats){
    for (j in i:seats){
      nodes[t,] = c(i,j-i,seats-j);
      t=t+1;
    }
  }  
  
  #allocate seats
  S  = matrix(apply(R, 1, function(x) alloc(letters[24:26], x, seats, step, threshold)), nrow=dots, ncol=3, byrow=TRUE)
  
  #Indexes for Voronoi and Allocation regions
  Uniform   = apply(R, 1, UniformNearest,   nodes = nodes, nnodes=nnodes)
  Manhattan = apply(R, 1, ManhattanNearest, nodes = nodes, nnodes=nnodes)
  Euclid    = apply(R, 1, EuclidNearest,    nodes = nodes, nnodes=nnodes)
  Allocated = apply(S, 1, AllocatedNode,    nodes = nodes, nnodes=nnodes)
  
  #node labels for seats
  label=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
  label=do.call("paste0",as.data.frame(label))
  label=rbind(matrix("",dots),as.matrix(label))
  
  #dates column
  time=matrix(NA,dots+nnodes)
  if (norm(dates)!=0) time = rbind(time, dates)
  
  #attach node points                
  x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/seats))
  y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/seats))
  z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/seats))
  
  t=0;
  if (norm(votes) != 0){
    
    x <- c(x, votes[,1])
    y <- c(y, votes[,2])
    z <- c(z, votes[,3])
    t=dim(votes)[1]
    label = c(label,as.character(dates))
    
  }
  
  #assembling a data frame
  df = data.frame(
    x ,
    y ,
    z ,
    
    code          = c(rgb(S/seats), matrix(NA,nnodes+t)),
    
    Euclid        = c(Euclid, matrix(NA,nnodes+t)),
    codeEuclid    = c(rgb(nodes[Euclid[1:dots],]/seats), matrix(NA,nnodes+t)),
    
    Manhattan     = c(Manhattan, matrix(NA,nnodes+t)),
    codeManhattan = c(rgb(nodes[Manhattan[1:dots],]/seats), matrix(NA,nnodes+t)),
    
    Uniform       = c(Uniform, matrix(NA,nnodes+t)),
    codeUniform   = c(rgb(nodes[Uniform[1:dots],]/seats), matrix(NA,nnodes+t)),
    
    Allocated     = c(Allocated, matrix(NA,nnodes+t)),
    Malapportionment = c((Allocated!=Euclid)[1:dots], matrix(NA,nnodes+t)),
    
    label,
    
    time
  )
  
  return(df);
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
alloc(letters[1:3], votes, 10, 1, .05) 


#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=20000
threshold=.05

#Generate random electoral results

#Ternary simulation. Projects points in[0,1]^3 to the sum(x)=1 hyperplane
#This creates way more points in the center than in the extrems
R = matrix(runif(3*dots), nrow=dots, ncol=3)
R = R/rowSums(R)

#Map Cartesian to Ternary to produce a homegeneous simulation
Rc = matrix(runif(2*dots), nrow=dots, ncol=2)
R  = matrix(apply(Rc,1, CartesianToTernary), nrow=dots, ncol=3, byrow=TRUE)

#Assemble the data frame
df = AssembleData(R, seats, step, threshold)


#Allocation
ggtern(data=df,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Allocation")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")

#Voronoi
ggtern(data=df,aes(x,y,z,color=codeManhattan, alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  geom_text(aes(label=label),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Voronoi")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")

#equivalence between Voronoi regions w/ different distances and between those and electoral regions
sum(df$Euclid[1:dots]==df$Manhattan[1:dots])/dots
sum(df$Uniform[1:dots]==df$Manhattan[1:dots])/dots
sum(df$Uniform[1:dots]==df$Euclid[1:dots])/dots
#points allocated in their corresponding Voronoi region
sum(df$Malapportionment[1:dots])/dots
#size of the regions # 1 ~= dots/nnodes
nnodes=(seats+1)*(seats+2)/2;
RegionSize = table(df$Allocated[1:dots])/(dots/nnodes)
plot(RegionSize)

#malapportionment
ggtern(data=df,aes(x,y,z,color=!Malapportionment, alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Malapportionment")+
  scale_colour_grey(na.value = "black")




###############################################
#history of election results
#example
#seats=7; step=1;
nelect=10;
votes=matrix(0,nelect,3)
votes[,1]=c(355,437,282,393,408,479,468,356,287,141)
votes[,3]=c(258,  0,185,144,159, 90, 47, 73,114,140)
votes[,2]=c(  0,  0, 82,  0,  0,  0,  0,  0,  0,139)
votes <- votes/rowSums(votes)
dates=as.matrix(seq(from=1979, to=2015, by=4))

dfvotes = AssembleData(R, seats, step, threshold, votes, dates)

#df2 = df[(dots+nnodes+1):(dots+nnodes+dim(dates)),]


ggtern(data=dfvotes,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  geom_text(aes(label=label), color="red", hjust=0.5,vjust=-0.6, size=3)+
  #geom_path(data=dfvotes,colour="blue", linetype=3, size=0.7)+
  labs(x="SocLib",y="SocCom",z="LibCon",title="Past Elections")+
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")
