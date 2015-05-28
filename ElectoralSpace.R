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

generateColors <- function (colorRGB, seats) {
  
  nnodes = (seats+1)*(seats+2)/2;
  nodes <- generateNodes(seats)
  #decimal codes for colors in nodes
  CodeDecRGB = floor(nodes*255.9/seats) #255.9 avoids seats -> 256 -> HEX #100 case
  #linear transformation to fixed extreme colors
  Code3 = floor(CodeDecRGB%*%matrix(colorRGB,3,3)/256)
  #hex codes for colors in nodes
  CodeRGB=cbind(matrix("#",nnodes),format(as.hexmode(Code3[,1]),width=2),format(as.hexmode(Code3[,2]),width=2),format(as.hexmode(Code3[,3]),width=2))
  CodeRGB=do.call("paste0",as.data.frame(CodeRGB))
  values=CodeRGB
  
  return(values);
  
}
  
AssembleData <- function (dots, seats, step, threshold, vectorOrderCode=as.matrix(3^(0:(seats-1))), votes=matrix(0,1), election=matrix(0,1), partial=seats) {
  
  #Generate random electoral results
  R <- generateDots(dots)
  
  #nodes
  nodes  = generateNodes(seats)
  nnodes = (seats+1)*(seats+2)/2;
  
  #allocate seats
  AllocStructure <- apply(R, 1, function(x) alloc(as.character(1:3), x, seats, step, threshold, partial=partial))
  #input partial=1 to compute partial sums or partial=seats to compute only seats
  
  #allocation
  AllocPartial=matrix("",dots,seats-partial+1)
  
  for (i in partial:seats){
    
    S = t(matrix(sapply(AllocStructure, function(x) x[[i-partial+1]]),nrow=3,ncol=dots))
    Code3 = floor(S*255.9/i) #255.9 avoids i (total seats) -> 256 -> HEX #100 case
    CodeRGB = cbind(matrix("#",dots),format(as.hexmode(Code3[,1]),width=2),format(as.hexmode(Code3[,2]),width=2),format(as.hexmode(Code3[,3]),width=2))
    CodeRGB = do.call("paste0", as.data.frame(CodeRGB))
    
    AllocPartial[,i-partial+1] = CodeRGB
    
  }
  
  AllocCode=AllocPartial[,seats-partial+1]
  
  Code3 = floor(S*255.9/seats) #255.9 avoids 7 -> 256 -> HEX #100 case
  CodeRGB = cbind(matrix("#",dots),format(as.hexmode(Code3[,1]),width=2),format(as.hexmode(Code3[,2]),width=2),format(as.hexmode(Code3[,3]),width=2))
  CodeRGB = do.call("paste0", as.data.frame(CodeRGB))
  
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
    
    Sx            = c(S[,1], matrix(NA,nnodes+el)),
    Sy            = c(S[,2], matrix(NA,nnodes+el)),
    Sz            = c(S[,3], matrix(NA,nnodes+el)),
    
    Allocated     = c(Allocated, matrix(NA,nnodes+el)),
    
    AllocCode     =  c(AllocCode, matrix(NA,nnodes+el)),
    
    CodeRGB         = c(CodeRGB, matrix(NA,nnodes+el)),
    
    AllocOrderCode=  c(AllocOrderCode, matrix(NA,nnodes+el)),
    
    Euclid        = c(Euclid, matrix(NA,nnodes+el)),
    Manhattan     = c(Manhattan, matrix(NA,nnodes+el)),
    Uniform       = c(Uniform, matrix(NA,nnodes+el)),
    
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
seats=7;
step=1; #(2 Sainte-Laguë 1 D'Hondt)
dots=50000
threshold=.05

#Assemble the data frame
#df = AssembleData(dots, seats, step, threshold)
df = AssembleData(dots, seats, step, threshold, partial=TRUE)

df2 =AssembleData(dots, seats, step=2,threshold, partial=2)

#alternative encoding for AllocOrder:
#default encoding:
#vectorOrderCoder = as.matrix(3^(0:(seats-1)))
#reverse encoding:
#vectorOrderCoderev = as.matrix(rev(3^(0:(seats-1))))
#df = AssembleData(dots, seats, step, threshold, vectorOrderCode=vectorOrderCoderev)

#colors for palettes
colorRGB1 <- c(242,74,87, 124,218,198, 91,168,246) #c("#f24a57","#7cdac6","#5ba8f6")
colorRGB2 <- c(226,10,23, 100,207,151, 0,112,184) #c("#e20a17","#64cf97","#0070b8")
colorRGB3 <- c(242,74,87, 124,218,198, 2,114,195)
colorRGB4 <- c(256,0,0, 0,256,0, 0,0,256)
colorRGB5  <- c(224,0,0, 96,224,176, 96,128,224)
colorRGB5  <- c(240,0,0, 96,224,176, 96,128,224) ## <- this one

colorRGB0 <- colorRGB5 #choose color palette

#Allocation

a1 <- ggtern(data=df,aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_text(aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="D'Hondt Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,seats), guide=FALSE, na.value="khaki2")

a2 <- ggtern(data=df2,aes(x,y,z,color=as.factor(Allocated)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_text(aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,seats), guide=FALSE, na.value="khaki2")

ggtern.multi(a1, a2, cols=2)

#Voronoi

a3 <- ggtern(data=df,aes(x,y,z,color=as.factor(Manhattan)))+
  theme_bw()+
  geom_point(alpha=1)+
  geom_text(aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="Voronoi Allocation")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB0,seats), guide=FALSE, na.value="khaki2")

ggtern.multi(a2, a3, cols=2)

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
m1 <- ggtern(data=df,aes(x,y,z,color=Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="D'Hondt Malapportionment")+
  #scale_colour_grey(na.value = "black", guide = FALSE)
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

m2 <- ggtern(data=df2,aes(x,y,z,color=Malapportionment)) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë Malapportionment")+
  scale_colour_brewer(palette = "YlGnBu", na.value = "grey60", guide = FALSE)

ggtern.multi(m1, m2, cols=2)

#ordering subregions
ggtern(data=df,aes(x,y,z,color=as.factor(AllocOrderCode))) +
  theme_rgbw() +
  geom_point(alpha=0.8) +
  labs(x="X",y="Y",z="Z",title="Allocation ordering regions")+
  scale_colour_grey(start = 0.1, end = 1, na.value = "black", guide = FALSE)

#Partial Allocations

p1 <- ggtern(data=df,aes(x,y,z,color=as.factor(df[,17+2])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 2 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,2), guide=FALSE)

p2 <- ggtern(data=df,aes(x,y,z,color=as.factor(df[,17+3])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 3 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,3), guide=FALSE)

p3 <- ggtern(data=df,aes(x,y,z,color=as.factor(df[,17+4])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="D'Hondt, 4 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,4), guide=FALSE)

p4 <- ggtern(data=df2,aes(x,y,z,color=as.factor(df2[,17+1])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 2 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,2), guide=FALSE)

p5 <- ggtern(data=df2,aes(x,y,z,color=as.factor(df2[,17+2])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 3 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,3), guide=FALSE)

p6 <- ggtern(data=df2,aes(x,y,z,color=as.factor(df2[,17+3])))+
  theme_bw()+
  geom_point(alpha=1)+
  labs(x="X",y="Y",z="Z",title="Sainte-Laguë, 4 seats")+
  scale_colour_manual(values=generateColors(colorRGB0,4), guide=FALSE)

plist=list(p1,p4,p2,p5,p3,p6)
ggtern.multi(plotlist=plist, cols=3)

#history of election results
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

#votes[,1]=c(355,437,282,393,408,479,468,356,287,234)
#votes[,3]=c(258,  0,185,144,159, 90, 47, 73,114, 75)
#votes[,2]=c(  0,  0, 82,  0,  0,  0,  0,  0,  0,138)

votes <- votes/rowSums(votes)
nnodes = (seats+1)*(seats+2)/2;

dfvotes = AssembleData(dots, seats, step, threshold, votes=votes, election=election)

dfSpline <- generateSpline(dfvotes, dots, seats, nelect)

ggtern(data=dfvotes,aes(x,y,z,color=as.factor(Allocated)))+ #CodeRGB also works
  #theme_rgbw()+
  theme_bw()+
  geom_point(alpha=1)+
  geom_text(data=dfvotes[(dots+nnodes+1):(dots+nnodes+nelect),],aes(label=label), color="slateblue4", hjust=1.2, vjust=0.3, size=3, angle=90)+ 
  geom_text(data=dfvotes[(dots+1):(dots+nnodes),],aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  geom_path(data=dfSpline,colour="darkgreen", linetype=1, size=1)+ 
  #geom_path(data=dfvotes[(dots+nnodes+1):(dots+nnodes+nelect),],colour="brown3", linetype=1, size=0.7)+ 
  labs(x="SocLib",y="SocCom",z="LibCon",title="Past Elections")+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  scale_colour_manual(values=generateColors(colorRGB1,seats), guide=FALSE, na.value="khaki2")


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


