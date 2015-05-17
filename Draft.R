#Install ggtern from an online CRAN repository
install.packages("ggtern")
#Load the ggtern library
library(ggtern)
#install.packages("ggplot2")
#library(ggplot2)
install.packages("e1071")
library(e1071)

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

#Allocation example (step=2 Sainte-LaguÃ«; step=1 D'Hondt)
votes <- sample(1:1000, 5) 
votes 
alloc(letters[1:5], votes, 10, 1) 

#presets
seats=5;
step=1; #(2 Sainte-LaguÃ« 1 D'Hondt)
dots=50000

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
  
  S[i,]= alloc(letters[24:26], R[i,], seats, step);
  
  
  #computing the nearest nodes to each point
  #Euclid[i]=which.min(sqrt(rowSums((nodes/seats-R[i,])^2)))
  Euclid[i]=which.min(apply((nodes/seats - rep(1,nnodes) %*% t(R[i,]))^2,1,sum))
  #sqrt not needed bc/ it's a monotone function
  Manhattan[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,sum))
  #Manhattan[i]=which.min(rowSums(abs(nodes/seats-R[i,])))
  Uniform[i]=which.min(apply(abs(nodes/seats - rep(1,nnodes) %*% t(R[i,])),1,max))
  Nodes[i]=which.min(apply(abs(nodes - rep(1,nnodes) %*% t(S[i,])),1,max))
}

#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(nodes[,1]/seats))
y <-rbind(as.matrix(R[,2]),as.matrix(nodes[,2]/seats))
z <-rbind(as.matrix(R[,3]),as.matrix(nodes[,3]/seats))

#assembling a data frame
df = data.frame(x ,
                y ,
                z ,
                code= c(rgb(S/seats),matrix("#000000",nnodes)),
                codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nnodes)),
                codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nnodes)),
                codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nnodes)),
                codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes)) )

#main plot
ggtern(data=df,aes(x,y,z,color=code, size = code=="#000000",alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation")+
  #scale_colour_hue()
  scale_colour_hue(l=25,c=180, h.start=0)

#ploting only the nodes
labeltext=matrix(NA,nnodes)
labeltext=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
#labeltext=cbind(as.character(x),matrix("-",nnodes),as.character(y),matrix("-",nnodes),as.character(z)),
labeltext=do.call("paste0",as.data.frame(labeltext))
df2 = data.frame(x=nodes[,1] ,
                 y=nodes[,2] ,
                 z=nodes[,3] ,
                 code= rgb(nodes/seats) )

ggtern(data=df2,aes(x,y,z)) +
  #theme_rgbw() +
  geom_point(size=2) +
  geom_text(aes(label=paste0(labeltext)),hjust=0.5,vjust=-0.6, size=3)+
  labs(x="X",y="Y",z="Z",title="Allocation nodes")

#ploting Voronoi regions:

ggtern(data=df,aes(x,y,z,color=codeEuclid, size = codeEuclid=="#000000", alpha=0.8)) +
  #ggtern(data=df,aes(x,y,z,color=codeManhattan, size = codeManhattan=="#000000", alpha=0.8)) +
  #ggtern(data=df,aes(x,y,z,color=codeUniform, size = codeUniform=="#000000", alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Euclid Voronoi Regions")+
  #scale_colour_brewer(palette="Paired")
  scale_colour_hue(l=25,c=180, h.start=0)#palette="Paired")

#equivalence between Voronoi regions w/ different distances and between those and electoral regions
sum(Euclid==Manhattan)/dots
sum(Uniform==Manhattan)/dots
sum(Uniform==Euclid)/dots
#points allocated in their corresponding Voronoi region
sum(Nodes==Euclid)/dots
#size of the regions # 1 ~= dots/nnodes
sapply(1:nnodes, function(x) sum(Nodes==x)*nnodes/dots)


#malapportionment
ggtern(data=df,aes(x,y,z,color=codeVoronoi, size = codeEuclid=="#000000", alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Malapportionment")#+
#scale_colour_brewer(palette="Paired")
#scale_colour_hue(l=25,c=180, h.start=0)#palette="Paired")


#scale_colour_hue(l=60,c=180, h.start=0)#palette="Paired") 
#scale_colour_hue(l=25,c=180, h.start=0)#palette="Paired")

################################3
#election result plots

#random results
votes <- sample(1:10000, 3) 
votes <- votes/sum(votes)
rS=alloc(letters[1:3], votes, 10, 1) 
rcode=matrix(1,nnodes)
#ploting only the nodes
df2 = data.frame(x=c(nodes[,1],votes[1]) ,
                 y=c(nodes[,2],votes[2]) ,
                 z=c(nodes[,3],votes[3]) , 
                 rcode=c(rcode,2) )


ggtern(data=df2,aes(x,y,z, color=rcode, alpha=0.9)) +
  #theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation nodes")


#example
#seats=7; step=1;
nelect=10;
votes=matrix(0,3,nelect)
votes[1,]=c(355,437,282,393,408,479,468,356,287,141)
votes[3,]=c(258,  0,185,144,159, 90, 47, 73,114,140)
votes[2,]=c(  0,  0, 82,  0,  0,  0,  0,  0,  0,139)
votes <- votes/sum(votes)

#ploting only the nodes
df3 = data.frame(x=votes[1,] ,
                 y=votes[2,] ,
                 z=votes[3,] ,
                 t=seq(from=1979, to=2015, by=4)
)
df3b = data.frame(x=votes[1,1:9] ,
                  y=votes[2,1:9] ,
                  z=votes[3,1:9]
)

ggtern(data=df3,aes(x,y,z, size=2, alpha=0.9)) +
  #theme_rgbw() +
  geom_point() +
  geom_path(data=df3b,colour="blue", linetype=3, size=0.7)+
  labs(x="SocLib",y="SocCom",z="LibCon",title="Allocation nodes")


##########
#results plus allocation plot #### plot a smooth curve through the results
#### not completed  # as a functiom of time


#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(votes[1,]))
y <-rbind(as.matrix(R[,2]),as.matrix(votes[2,]))
z <-rbind(as.matrix(R[,3]),as.matrix(votes[3,]))
#assembling a data frame
df4 = data.frame(x,
                 y,
                 z,
                 #code= c(S[,3],matrix(NA,nelect)),
                 code= c(rgb(S/seats),matrix(NA,nelect)),
                 codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nelect)),
                 codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nelect)),
                 codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nelect)),
                 codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nelect)),
                 labeltext=c(matrix("",dots),seq(from=1979, to=2015, by=4))
)


ggtern(data=df4,aes(x,y,z,colour=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point(size=3) +
  geom_text(aes(label=as.character(labeltext)),hjust=-0.1,just=0, size=4)+
  geom_path(data=df3,colour="blue", linetype=3, size=0.7)+
  labs(x="SD",y="SOC",z="CL",title="Past Elections")+
  #scale_colour_brewer()
  #scale_colour_hue()
  #scale_colour_hue(l=70,c=130, h.start=50)
  scale_colour_grey(start = 0.4, end = 1, na.value = "red")

#geom_text(aes(label=ifelse(code==NA,"a",'')),hjust=0,just=0)+

############################
#attaching seats numbers to the diagram

escala=matrix(0,seats+1,4)
escala[,1]=seq(from=1/(2*(seats+1)), to=1, by=1/(seats+1))
escala[,2]=seq(from=(1-1/(2*(seats+1)))/2, to=0, by=-1/(2*(seats+1)))
escala[,3]=seq(from=(1-1/(2*(seats+1)))/2, to=0, by=-1/(2*(seats+1)))
escala[,4]=0:seats
rowSums(escala)
escala[,c(2,1,3,4)]
escala[,c(2,3,1,4)]
escalas = rbind(escala, escala[,c(2,1,3,4)],escala[,c(2,3,1,4)])

#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(votes[1,]),as.matrix(escalas[,1]))
y <-rbind(as.matrix(R[,2]),as.matrix(votes[2,]),as.matrix(escalas[,2]))
z <-rbind(as.matrix(R[,3]),as.matrix(votes[3,]),as.matrix(escalas[,3]))
t <-rbind(matrix(0,dots),as.matrix(seq(from=1979, to=2015, by=4)),as.matrix(escalas[,4]))

#assembling a data frame
df4b = data.frame(x,
                  y,
                  z,
                  t,
                  #code= c(S[,3],matrix(NA,nelect)),
                  code=c(rgb(S/seats),matrix(NA,nelect),matrix("#000000",3*(seats+1))),
                  #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nelect)),
                  #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nelect)),
                  #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nelect)),
                  #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nelect)),
                  labeltext=c(matrix("",dots),seq(from=1979, to=2015, by=4), escalas[,4] ) 
)

ggtern(data=df4b,aes(x,y,z,colour=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point(size=3) +
  geom_text(aes(label=as.character(labeltext)),hjust=-0.2,vjust=-0.2, size=4)+
  #geom_path(data=df3,colour="blue", linetype=3, size=0.7)+
  geom_smooth(data=df3, method="loess", formula=y~x,se=F,limitarea=F,fullrange=T,
              color="blue",size=1,linetype=5)+
  labs(x="SD",y="SOC",z="CL",title="Past Elections")+
  #scale_colour_brewer()
  #scale_colour_hue()
  #scale_colour_hue(l=70,c=130, h.start=50)
  scale_colour_grey(start = 0.4, end = 1, na.value = "red")


#plot(df3$t,df3$x, type='b', col='red')
#xspline(df3$t,df3$x, shape=1)


##########classification  ###not completed

dat=data.frame(code=df5$code[1:dots],
               x=df5$x[1:dots],
               y=df5$y[1:dots],
               z=df5$z[1:dots]
)

fit <- svm(formula = code~., data = dat, kernel = "radial", cost=10, gamma=1)
print(fit)
plot(fit)

############################
#attaching seat levels to the diagram

escala=matrix(0,seats+1,4)
escala[,1]=seq(from=1/(2*(seats+1)), to=1, by=1/(seats+1))
escala[,2]=seq(from=(1-1/(2*(seats+1)))/2, to=0, by=-1/(2*(seats+1)))
escala[,3]=seq(from=(1-1/(2*(seats+1)))/2, to=0, by=-1/(2*(seats+1)))
escala[,4]=0:seats
rowSums(escala)
escala[,c(2,1,3,4)]
escala[,c(2,3,1,4)]
escalas = rbind(escala, escala[,c(2,1,3,4)],escala[,c(2,3,1,4)])

#attaching node points                
x <-rbind(as.matrix(R[,1]),as.matrix(escalas[,1]))
y <-rbind(as.matrix(R[,2]),as.matrix(escalas[,2]))
z <-rbind(as.matrix(R[,3]),as.matrix(escalas[,3]))
t <-rbind(matrix(0,dots),as.matrix(escalas[,4]))

#assembling a data frame
df5 = data.frame(x,
                 y,
                 z,
                 t,
                 #code= c(S[,3],matrix(NA,nelect)),
                 code=c(rgb(S[,c(2,1,3)]/seats),matrix(NA,3*(seats+1))),
                 #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nelect)),
                 #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nelect)),
                 #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nelect)),
                 #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nelect)),
                 labeltext=c(matrix("",dots), escalas[,4] ) 
)

ggtern(data=df5,aes(x,y,z,colour=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point(size=2) +
  geom_text(aes(label=as.character(labeltext)),hjust=-0.2,vjust=-0.2, size=4)+
  labs(x="X",y="Y",z="Z",title="Allocating 7 seats to 3 parties")+
  #scale_colour_brewer()
  #scale_colour_hue()
  #scale_colour_hue(l=70,c=130, h.start=50)
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")


###################################
###################################
#plot regions and nodes #definitive diagram
label1=matrix(NA,dots)
label2=label1
labeltext=label2
label1=cbind(as.character(nodes[,1]),matrix("-",nnodes),as.character(nodes[,2]),matrix("-",nnodes),as.character(nodes[,3]))
label2=do.call("paste0",as.data.frame(label1))
labeltext=rbind(matrix(NA,dots),as.matrix(label2))
df6=data.frame(x=rbind(as.matrix(R[,1]),as.matrix(nodes[,1])) ,
               y=rbind(as.matrix(R[,2]),as.matrix(nodes[,2])) ,
               z=rbind(as.matrix(R[,3]),as.matrix(nodes[,3])) ,
               
               #labeltext=cbind(as.character(x),matrix("-",nnodes),as.character(y),matrix("-",nnodes),as.character(z)),
               
               labeltext,
               code= c(rgb(S[,c(2,1,3)]/seats),matrix(NA,nnodes))  )

ggtern(data=df6,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  geom_text(aes(label=labeltext),hjust=0.5,vjust=-0.6, size=5)+
  labs(x="X",y="Y",z="Z",title="Allocating 5 seats to 3 parties")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")

###################################
#alternative series of divisors

seats=7
div=matrix(1,seats)
#define a series of divisors
for (i in 3:seats){
  div[i]=div[i-1]+div[i-2]
}
#
step=sqrt(2)
for (i in 2:seats){
  div[i]=div[i-1]*step
}

alloc2 <- function(parties, votes, seats){ 
  qtable <- data.frame( 
    
    parties = rep(parties, each = seats), 
    quotients  = as.vector(sapply(votes, function(x) x/div)),
    votesrep = rep(votes, each = seats)
  ) 
  qtable <- qtable$parties[order(-qtable$quotients, -qtable$votesrep)] [1:seats] 
  table(qtable) 
}
votes <- sample(1:1000, 5) 
votes 
alloc2(letters[1:5], votes, seats) 

#################################################
#################################################
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
  #do.call("paste0",as.data.frame(as.matrix(qtable)))
  
}


#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 5) 
votes
#votes/sum(votes)
allocOrder(letters[1:5], votes, 10, 1, .05) #let's name parties with integers


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



##################################################
#ordering regions

#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=30000
threshold=.05
vectorcode=matrix(1,seats)
#for (i in 2:seats)
#  vectorcode[i]=vectorcode[i-1]*3

#vectorcode
vectorcode = as.matrix(3^(0:(seats-1)))
vectorcoderev = as.matrix(rev(3^(0:(seats-1))))

#defining a function for the seat allocation order with an entry threshold
allocOrderCode <- function(parties, votes, seats, step, threshold){ 
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
  #do.call("paste0",as.data.frame(as.matrix(qtable)))
  t(as.matrix(qtable)-1)%*%vectorcode
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
#votes/sum(votes)
allocOrderCode(1:3, votes, seats, 1, .05) #let's name parties with integers


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
OrderCode=matrix(0,dots)
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
  
  OrderCode[i]= allocOrderCode(1:3, R[i,], seats, step, threshold);
  
  
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
dfOrder = data.frame(x ,
                     y ,
                     z ,
                     code= c(rgb((OrderCode)/(3^seats),0,0),matrix(NA,nnodes))
                     #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nnodes)),
                     #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nnodes)),
                     #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nnodes)),
                     #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes))
)

#main plot
ggtern(data=dfOrder,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black")
  scale_colour_grey(start = 0, end = 1, na.value = "black")

###### reverse encoding ###

#defining a function for the seat allocation order with an entry threshold
allocOrderCodeRev <- function(parties, votes, seats, step, threshold){ 
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
  #do.call("paste0",as.data.frame(as.matrix(qtable)))
  t(as.matrix(qtable)-1)%*%vectorcoderev
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
#votes/sum(votes)
allocOrderCodeRev(1:3, votes, seats, 1, .05) #let's name parties with integers


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
OrderCodeRev=matrix(0,dots)
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
  
  OrderCodeRev[i]= allocOrderCodeRev(1:3, R[i,], seats, step, threshold);
  
  
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
dfOrderRev = data.frame(x ,
                        y ,
                        z ,
                        code= c(rgb((OrderCodeRev)/(3^seats),0,0),matrix(NA,nnodes))
                        #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nnodes)),
                        #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nnodes)),
                        #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nnodes)),
                        #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes))
)

#main plot
ggtern(data=dfOrderRev,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black")
  scale_colour_grey(start = 0.1, end = 1, na.value = "black")
#this encoding is more coherent but we can view the regions better with the previous one


#####################################################
##let's try a different encoding for ordering regions

#presets
seats=5;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=30000
threshold=.05
vectorcode2=matrix(1,seats)
#for (i in 2:seats)
#  vectorcode2[i]=vectorcode2[i-1]*2

#vectorcode2

vectorcode2=as.matrix(2^(0:(seats-1)))

#defining a function for the seat allocation order with an entry threshold
allocOrderCode2 <- function(parties, votes, seats, step, threshold){ 
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
  #do.call("paste0",as.data.frame(as.matrix(qtable)))
  #t(as.matrix(qtable)-1)%*%vectorcode
  cbind(t(as.matrix(qtable)==1)%*%vectorcode2,t(as.matrix(qtable)==2)%*%vectorcode2,t(as.matrix(qtable)==3)%*%vectorcode2)
  
}

#Allocation example (step=2 Sainte-Laguë; step=1 D'Hondt)
votes <- sample(1:1000, 3) 
votes
#votes/sum(votes)
allocOrderCode2(1:3, votes, seats, 1, .05) #let's name parties with integers


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
OrderCode2=matrix(0,dots,3)
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
  
  OrderCode2[i,]= allocOrderCode2(1:3, R[i,], seats, step, threshold);
  
  
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
dfOrder2 = data.frame(x ,
                      y ,
                      z ,
                      code= c(rgb(OrderCode2/(2^seats)),matrix(NA,nnodes))
                      #codeEuclid=c(rgb(nodes[Euclid,]/seats),matrix("#000000",nnodes)),
                      #codeManhattan=c(rgb(nodes[Manhattan,]/seats),matrix("#000000",nnodes)),
                      #codeUniform=c(rgb(nodes[Uniform,]/seats),matrix("#000000",nnodes)),
                      #codeVoronoi=c(Nodes==Euclid,matrix(FALSE,nnodes))
)

#main plot
ggtern(data=dfOrder2,aes(x,y,z,color=code,alpha=0.8)) +
  theme_rgbw() +
  geom_point() +
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  #scale_colour_hue()
  #scale_colour_hue(l=25,c=180, h.start=0)
  scale_colour_grey(start = 0.4, end = 1, na.value = "black")
#although this encoding is more formal I like the first diagram better
