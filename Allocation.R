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

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 5) 
votes 
alloc(letters[1:5], votes, 10, 1) 


#defining a function for the seat allocation with an entry threshold
allocT <- function(parties, votes, seats, step, threshold){ 
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


#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 9) 
votes
#votes/sum(votes)
allocT(letters[1:9], votes, 100, 1, .05) 



#defining a function for the seat allocation order with an entry threshold
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


#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 5) 
votes
#votes/sum(votes)
allocOrder(letters[1:5], votes, 10, 1, .05) 


#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)

#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=30000
threshold=.2

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
S=matrix(0,dots,3)
R=matrix(0,dots,3)
x=matrix(0,dots)
y=x
z=x
Euclid=matrix(0,dots)
Manhattan=matrix(0,dots)
Uniform=matrix(0,dots)
Nodes=matrix(0,dots)
R2=matrix(0,dots,2);
#R2=runif(dosts*2)
#dim(R2)<-c(dots,2)
for (i in 1:dots) {
  
  #first simulation. This creates way more points in the center
  #than in the extrems
  #this projects points in[0,1]^3 to the sum(x)=1 hyperplane
  #R[i,] = runif(3);
  #R[i,]=R[i,]/sum(R[i,]);
  ###########################
  
  #second simulation
  #R2[i,]=runif(2);
  #if (sum(R2[i,])>1) R2[i,]=matrix(1,2)-R2[i,]
  #R[i,] = cbind(R[i,1],R[i,2],1-R[i,1]-R[i,2]);
  #R[i,]=R[i,]/sum(R[i,]);
  ###########################
  
  #homogeneous simulation
  #send an subset of [0,1]x[0,1] which is an equilateral triangle to the ternary diagram 
  R2[i,]=runif(2);
  while ((R2[i,1]<=1/2 & R2[i,2]>2*sin(pi/3)*R2[i,1]) | (R2[i,1]>1/2 & R2[i,2]>2*sin(pi/3)*(1-R2[i,1])))
    R2[i,]=runif(2);
  
  
  R[i,] = cbind(1-R2[i,1]-R2[i,2]/(2*sin(pi/3)),R2[i,2]/sin(pi/3),R2[i,1]-R2[i,2]/(2*sin(pi/3)));
  #R[i,]=R[i,]/sum(R[i,]);
  ## here ends the third simulation
  ############################
  
  S[i,]= allocT(letters[24:26], R[i,], seats, step, threshold);
  
  
  #computing the nearest nodes to each point
 
  #Euclid[i]=which.min(apply((nodes/seats - rep(1,nnodes) %*% t(R[i,]))^2,1,sum))
  #sqrt not needed bc/ it's a monotone function
  #Manhattan[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,sum))
  #Uniform[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,max))
  #Nodes[i]=which.min(apply(abs(nodes - rep(1,nnodes) %*% t(S[i,])),1,max))
}

#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/seats))
y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/seats))
z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/seats))

#assembling a data frame
dfT = data.frame(x ,
                y ,
                z ,
                code= c(rgb(S[,c(2,1,3)]/seats),matrix(NA,nnodes))
                #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nnodes)),
                #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nnodes)),
                #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nnodes)),
                #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes))
                )

#main plot
ggtern(data=dfT,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation with threshold")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")
