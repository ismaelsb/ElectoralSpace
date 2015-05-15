#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)

#defining a function for the seat allocation
alloc <- function(parties, votes, seats, step){ 
  qtable <- data.frame( 
    parties = rep(parties, each = seats), 
    quotients  = as.vector(sapply(votes, function(x) x / 
                                    seq(from=1, to=1+step*(seats-1), by=step) )),
    votesrep = rep(votes, each = seats)
  ) 
  qtable <- qtable$parties[order(-qtable$quotients, -qtable$votesrep)] [1:seats] 
  table(qtable)
  
}

CartesianToTernary <- function (x) {
  
  while ((x[1]<=1/2 & x[2]>2*sin(pi/3)*x[1]) | (x[1]>1/2 & x[2]>2*sin(pi/3)*(1-x[1])))
    x=runif(2);
  
  v = cbind(1-x[1]-x[2]/(2*sin(pi/3)),x[2]/sin(pi/3),x[1]-x[2]/(2*sin(pi/3)));
  
  return(v);
  
}

uniformNode <- function (x) {
  
  Uniform=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(x)),1,max))
  
  return(Uniform);
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 5) 
votes 
alloc(letters[1:5], votes, 10, 1) 

#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=10000

#define nodes
t=1;
nnodes=(seats+1)*(seats+2)/2;
nodes= matrix(0,nnodes,3);
for (i in 0:seats){
  for (j in i:seats){
    nodes[t,] = c(i,j-i,seats-j);
    t=t+1;
  }
}

#simulation
R=matrix(runif(3*dots), nrow=dots, ncol=3)
S=matrix(0,dots,3)

x=matrix(0,dots)
y=x
z=x
Euclid=matrix(0,dots)
Manhattan=matrix(0,dots)
Uniform=matrix(0,dots)
Nodes=matrix(0,dots)
Rc=matrix(0,dots,2);


Rc=matrix(runif(2*dots), nrow=dots, ncol=2)
R = matrix(apply(Rc,1, CartesianToTernary), nrow=dots, ncol=3, byrow=TRUE)
Uniform=apply(R, 1, uniformNode)
S = matrix(apply(R, 1, function(x) alloc(letters[24:26], x, seats, step)), nrow=dots, ncol=3, byrow=TRUE)
 
  
  #Euclid[i]=which.min(apply((nodes/seats - rep(1,nnodes) %*% t(R[i,]))^2,1,sum))
  #Manhattan[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,sum))
  #Uniform[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,max))
  #Nodes[i]=which.min(apply(abs(nodes - rep(1,nnodes) %*% t(S[i,])),1,max))


#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/seats))
y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/seats))
z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/seats))

#assembling a data frame
df = data.frame(x ,
                y ,
                z ,
                code= c(rgb(S/seats),matrix(NA,nnodes))
                #, codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix(NA,nnodes)),
                #, codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix(NA,nnodes)),
                #, codeUniform=c(rgb(nodes[Uniform,]/seats),matrix(NA,nnodes))
                #, codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes)) 
)

#main plot
ggtern(data=df,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")
