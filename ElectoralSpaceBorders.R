#install.packages("ggtern")
library(ggtern)


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

generateBorderLattice <- function ( seats, step=1 ) {
  
  divisor <- seq(from=1, to=1+step*(seats+2-1), by=step)
  divisor <- c(0,divisor)
  
  counter <- 1
  borderpoint <- matrix(0,nrow=(seats+3)^2,ncol=4)
  # (seats+2)*(seats+3)/2 + (seats+3)*(seats+4)/2 = (seats+3)^2
  # (corners are counted and registered twice)
  
  #every sum seats+1 border points -> (seats+2)*(seats+3)/2
  #surplus 1 points
  for (i in 0:(seats+1)){
    
    for (j in 0:(seats+1-i)){
      
      d1 <- divisor[i+1] # +1 because vector element divisor[1]=0
      d2 <- divisor[j+1]
      d3 <- divisor[seats+1-i-j+1]
      
      borderpoint[counter,] <- matrix(c(d1,d2,d3,1))
      counter <- counter+1
      
    }
    
  }
  
  #every sum seats+2 border points -> (seats+3)*(seats+4)/2
  #surplus 2 points
  for (i in 0:(seats+2)){
    
    for (j in 0:(seats+2-i)){
      
      d1 <- divisor[i+1]
      d2 <- divisor[j+1]
      d3 <- divisor[seats+2-i-j+1]
      
      borderpoint[counter,] <- matrix(c(d1,d2,d3,2))
      counter <- counter+1
      
    }
    
  }
  
  
  
  borderpoint <- as.data.frame(borderpoint)
  names(borderpoint) <- c('x','y','z','surplus')
  
  return(borderpoint)
  
}

generateBorderLines <- function ( seats, step=1 ) {
  
  dB <- generateBorderLattice(seats,step)
  Bl <- matrix(0,nrow=3*(seats+1)*(seats+2)/2+3*(seats+1)+3,ncol=6)
        #number of edges: 6 per hexagon (counted twice);
        #                 2 per external hexagon (counted once)
        #                 3 autoedges on the corners (needed)
  
  #edges form each surplus 1 point to the 3 adjacent suplus 2 points
  
  Bl[,1] <- rep(dB[dB$surplus==1,1],each=3)
  Bl[,2] <- rep(dB[dB$surplus==1,2],each=3)
  Bl[,3] <- rep(dB[dB$surplus==1,3],each=3)
  
  Id3 <- matrix(rep(c(1,0,0,0,1,0,0,0,1),(seats+2)*(seats+3)/2),ncol=3, byrow=T)
  
  # surplus 1 to adjacent surplus 2
  extToint <- (Bl[,1:3]!=0)*Id3  #remove external points growing to internal points
  Bl[,4:6] <- Bl[,1:3]+Id3+(step-1)*extToint
              # external points grow by 1 to inernal points
              # external points grow by step to external points
              # internal points always grow by step
  
  Bl <- as.data.frame(Bl)
  names(Bl) <- c('x','y','z', 'xend','yend','zend') #extrem points for an edge
  
  return(Bl)
  
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


seats=5
step=1

NodesData <- generateNodes(seats)

dfB <- generateBorderLattice(seats, step)

dfB[,1:3]/rowSums(dfB[,1:3]) #integer to racional

dfBl <- generateBorderLines(seats, step)

ggtern( data=dfB, aes(x,y,z))+ 
  geom_point()

bl1 <- ggtern( data=dfBl, aes(x,y,z,xend=xend,yend=yend,zend=zend))+ 
  geom_segment()

ggtern(data=dfB, aes(x,y,z))+
  geom_point()+
  geom_segment(data=dfBl,aes(x,y,z,xend=xend,yend=yend,zend=zend))

bpl1 <- ggtern(data=dfB,aes(x,y,z))+
  #theme_bw()+
  geom_point(alpha=1, color=dfB$surplus)+
  geom_segment(data=dfBl,aes(x,y,z,xend=xend,yend=yend,zend=zend))+
  geom_point(data=NodesData,aes(x,y,z),color="khaki2")+
  geom_text(data=NodesData,aes(label=label), color="grey30", hjust=0.5, vjust=-0.6, size=4)+ 
  labs(x="X",y="Y",z="Z",title="D'Hondt Allocation")#+
  #geom_polygon(aes(fill=, group=))+
  #scale_colour_grey(start = 0.4, end = 1, na.value = "black", guide = FALSE)
  #scale_colour_manual(values=generateColors(colorRGB0,seats), guide=FALSE, na.value="khaki2")


step=2

dfB <- generateBorderLattice(seats, step)

#dfB[,1:3]/rowSums(dfB[,1:3]) #integer to racional

dfBl <- generateBorderLines(seats, step)

bl2 <- ggtern( data=dfBl, aes(x,y,z,xend=xend,yend=yend,zend=zend))+ 
  geom_segment()


ggtern.multi(bl1, bl2, cols=2)
bpl1

bl <- NULL

for (i in 1:4){
  
  dfBl <- generateBorderLines(seats=i, step=1)
  bl[[i]] <- ggtern( data=dfBl, aes(x,y,z,xend=xend,yend=yend,zend=zend))+ 
    geom_segment()
  
}

for (i in 1:4){
  
  dfBl <- generateBorderLines(seats=i, step=2)
  bl[[i+4]] <- ggtern( data=dfBl, aes(x,y,z,xend=xend,yend=yend,zend=zend))+ 
    geom_segment()
  
}

ggtern.multi(bl[[1]],bl[[5]],bl[[2]],bl[[6]],bl[[3]],bl[[7]],bl[[4]],bl[[8]], cols=4)




# Effective thresholds

EffThresholdTable <- function (maxseats, step=1){
  
  dfBlist <- NULL
  dfBlistQ <- NULL
  
  for (i in 1:maxseats){
    
    dfBlist[[i]] <- generateBorderLattice(seats=i, step)            #integers
    dfBlistQ[[i]] <- dfBlist[[i]][,1:3]/rowSums(dfBlist[[i]][,1:3]) #rationals
    
  }
  
  EffNecThreshold <- NULL
  EffSufThreshold <- NULL
  
  for (i in 1:maxseats){
    
    nonzero1 <- dfBlist[[i]][,1]!=0 & dfBlist[[i]]$surplus==1
    nonzero2 <- dfBlist[[i]][,1]!=0 & dfBlist[[i]]$surplus==2
    
    EffNecThreshold[i] <- min(dfBlistQ[[i]][nonzero1,1])
    EffSufThreshold[i] <- min(dfBlistQ[[i]][nonzero2,1])
    
  }
  
  EffThreshold <- as.data.frame(cbind(1:maxseats,EffNecThreshold, EffSufThreshold))
  names(EffThreshold) <- c('seats','EffNecThreshold','EffSufThreshold')
  
  return(EffThreshold)
  
}


maxseats <- 32
step=1

EffThreshold <- EffThresholdTable(maxseats, step)

ggplot(EffThreshold, aes(seats, EffNecThreshold))+
  #geom_point(colour='paleturquoise4')+
  scale_x_continuous(breaks=seq(1,maxseats,2))+
  scale_y_continuous(breaks=seq(0,0.5,0.02))+
  geom_line(aes(y=EffNecThreshold), colour='black', alpha=1/3)+
  geom_line(aes(y=EffSufThreshold), colour='black', alpha=1/2)+
  geom_ribbon(aes(ymin=EffNecThreshold,ymax=EffSufThreshold),fill='paleturquoise4', alpha=.85)+
  labs(title='Necessary votes for representation.\n A lighter shade for a greater number of parties, from two in the upper curve.\n Upper curve also for sufficient votes, for any number of parties',
       x='Total number of seats', y='Representation thresholds')+
  #geom_line(aes(y=1/(seats+1)),colour='black', alpha=1/2)+ # suff for 3 parties, already painted
  #geom_line(aes(y=1/(seats+2)),colour='black', alpha=1/3)+  # nec for 3 parties, already painted
  geom_line(aes(y=1/(seats+3)),colour='black', alpha=1/4)+ #theoretical nec thresholds from 3 parties
  geom_line(aes(y=1/(seats+4)),colour='black', alpha=1/5)+
  geom_line(aes(y=1/(seats+5)),colour='black', alpha=1/6)+
  geom_line(aes(y=1/(seats+6)),colour='black', alpha=1/7)+
  geom_line(aes(y=1/(seats+7)),colour='black', alpha=1/8)+
  geom_line(aes(y=1/(seats+8)),colour='black', alpha=1/9)+
  geom_line(aes(y=1/(seats+9)),colour='black', alpha=1/10)+
  geom_line(aes(y=1/(seats+10)),colour='black', alpha=1/11)+
  geom_line(aes(y=1/(seats+11)),colour='black', alpha=1/12)+
  geom_line(aes(y=1/(seats+12)),colour='black', alpha=1/13)+
  geom_line(aes(y=1/(seats+13)),colour='black', alpha=1/14)+
  geom_line(aes(y=1/(seats+14)),colour='black', alpha=1/15)+
  #geom_ribbon(aes(ymin=1/(seats+2),ymax=1/(seats+1)),fill='paleturquoise4', alpha=.85)+
  geom_ribbon(aes(ymin=1/(seats+3),ymax=1/(seats+2)),fill='paleturquoise4', alpha=.85^2)+
  geom_ribbon(aes(ymin=1/(seats+4),ymax=1/(seats+3)),fill='paleturquoise4', alpha=.85^3)+
  geom_ribbon(aes(ymin=1/(seats+5),ymax=1/(seats+4)),fill='paleturquoise4', alpha=.85^4)+
  geom_ribbon(aes(ymin=1/(seats+6),ymax=1/(seats+5)),fill='paleturquoise4', alpha=.85^5)+
  geom_ribbon(aes(ymin=1/(seats+7),ymax=1/(seats+6)),fill='paleturquoise4', alpha=.85^6)+
  geom_ribbon(aes(ymin=1/(seats+8),ymax=1/(seats+7)),fill='paleturquoise4', alpha=.85^7)+
  geom_ribbon(aes(ymin=1/(seats+9),ymax=1/(seats+8)),fill='paleturquoise4', alpha=.85^8)+
  geom_ribbon(aes(ymin=1/(seats+10),ymax=1/(seats+9)),fill='paleturquoise4', alpha=.85^9)+
  geom_ribbon(aes(ymin=1/(seats+11),ymax=1/(seats+10)),fill='paleturquoise4', alpha=.85^10)+
  geom_ribbon(aes(ymin=1/(seats+12),ymax=1/(seats+11)),fill='paleturquoise4', alpha=.85^11)+
  geom_ribbon(aes(ymin=1/(seats+13),ymax=1/(seats+12)),fill='paleturquoise4', alpha=.85^12)+
  geom_ribbon(aes(ymin=1/(seats+14),ymax=1/(seats+13)),fill='paleturquoise4', alpha=.85^13)+
  geom_ribbon(aes(ymin=0,ymax=1/(seats+14)),fill='paleturquoise4', alpha=.85^14)+
  theme()



maxseats <- 32
step=2

EffThreshold <- EffThresholdTable(maxseats, step)

ggplot(EffThreshold, aes(seats, EffNecThreshold))+
  #geom_point(colour='paleturquoise4')+
  scale_x_continuous(breaks=seq(1,maxseats,2))+
  scale_y_continuous(breaks=seq(0,0.5,0.02))+
  geom_line(aes(y=EffNecThreshold), colour='black', alpha=1/2)+
  geom_line(aes(y=EffSufThreshold), colour='black', alpha=1/2)+
  geom_ribbon(aes(ymin=EffNecThreshold,ymax=EffSufThreshold),fill='paleturquoise4', alpha=.85)+
  labs(title='Sainte-Laguë necessary votes for representation.\n Upper curve for sufficient votes.',
       x='Total number of seats', y='Representation thresholds')+
  theme()




#3Dplots

generateBorder4Lattice

generateBorder4Lines


