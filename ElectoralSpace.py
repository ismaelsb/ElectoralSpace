# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches


def generateHexagonData (seats, step=1):
    
    divisor = np.arange(1, 1+step*(seats+2-1), step)
    divisor = np.append(0,divisor)
    
    NodesData = generateNodes(seats)
    
    hexagon = np.zeros(((seats+1)*(seats+2)//2, 6*3),dtype=int)
    
    for i in range((seats+1)*(seats+2)//2):
        
        n = NodesData.ix[i,'x':'z'] #node
        hexagon[i,] = np.tile(np.array([n[0],n[1],n[2]],dtype=int), 6) + np.array([1,0,0, 1,0,1, 0,0,1, 0,1,1, 0,1,0, 1,1,0])    
        hexagon[i,] = divisor[hexagon[i,]]
    
    
    dhex = pd.concat([pd.DataFrame(np.repeat(range((seats+1)*(seats+2)//2),6)),
                      pd.DataFrame(np.reshape(hexagon,(6*(seats+1)*(seats+2)/2,3))),
                      pd.DataFrame(np.repeat(np.array(NodesData.x, dtype=int),6)),
                      pd.DataFrame(np.repeat(np.array(NodesData.y, dtype=int),6)),
                      pd.DataFrame(np.repeat(np.array(NodesData.z, dtype=int),6))],axis=1)
    
    dhex.columns = ['hex','x','y','z','nx','ny','nz']
    
    return(dhex)
  
 
def generateNodes (seats):
    
    nnodes = (seats+1)*(seats+2)//2
    nodes  = np.zeros((nnodes,3))
    t = 0
    
    for i in range(seats+1):
        
        for j in range(i,seats+1):
            
            nodes[t,] = [i,j-i,seats-j]
            t=t+1
            
    label = pd.concat([pd.DataFrame(nodes[:,0]), #pd.DataFrame(['-']*nnodes),
                       pd.DataFrame(nodes[:,1]), #pd.DataFrame(['-']*nnodes),
                       pd.DataFrame(nodes[:,2])], axis=1)
    label.columns = ['x','y','z']
    label=label.x.map(int).map(str)+'-'+label.y.map(int).map(str)+'-'+label.z.map(int).map(str)
    nodes = pd.concat([pd.DataFrame(nodes),label],axis=1)
    nodes.columns = ['x','y','z','label']
    
    return(nodes)


def generateColors (colorRGB, seats):
  
    ColorDec3 = np.array(colorRGB)/255    
    nodes = generateNodes(seats)[["x","y","z"]]
    #decimal codes for colors in nodes
    Code3 = nodes/seats
    #linear transformation to fixed extreme colors
    CodeDecRGB = np.dot(Code3, np.array(ColorDec3).reshape((3,3)))
    
    return(CodeDecRGB)

 
def ternaryToCartesian (data):
    
    TtoC = np.array([[0,0],
                     [1/2,np.sqrt(3)/2],
                     [1,0]])
                     
    return(np.dot(data,TtoC))
 

def plotallocation (seats, step=1, colorRGB0=[242,74,87, 124,218,198, 2,114,195], alpha=0.75):
  
    NodesData = generateNodes(seats)
    dhex = generateHexagonData(seats, step)
    nnodes = (seats+1)*(seats+2)//2
  
    dpolygnorm = dhex
    dpolygnorm[['x','y','z']] = dhex[['x','y','z']].div(dhex[['x','y','z']].sum(axis=1), axis=0)
    dpolygcart = pd.DataFrame(ternaryToCartesian(dpolygnorm[['x','y','z']]))
    dpolygcart.columns = ['x1','x2']
    dpolygnorm = pd.concat([dpolygnorm,dpolygcart],axis=1)
    
    dnodesnorm = NodesData
    dnodesnorm[['x','y','z']] = NodesData[['x','y','z']].div(NodesData[['x','y','z']].sum(axis=1), axis=0)
    dnodescart = pd.DataFrame(ternaryToCartesian(dnodesnorm[['x','y','z']]))   
    dnodescart.columns = ['x1','x2']
    dnodesnorm = pd.concat([dnodesnorm,dnodescart],axis=1)    
    
    CodeDecRGB = generateColors (colorRGB0, seats)
    
    fig = plt.figure(figsize=(7,7), dpi=160)
    ax = fig.add_subplot(111)

    #hexagon path instructions
    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    for i in range(nnodes):
        
        #hexagon vertices    
        verts = dpolygnorm.loc[dpolygnorm['hex'] == i][['x1','x2']] #select i-th polygon
            
        verts.set_index([[0,1,2,3,4,5]], inplace=True) #rename indexes
        verts.loc[6] = verts.loc[0]  #close the polygon 
        
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor=CodeDecRGB[i,], edgecolor='#404040', lw=.5, alpha=.95, zorder=0)
        ax.add_patch(patch)
        ax.annotate(dnodesnorm.label[i],xy=dpolygnorm.loc[dpolygnorm['hex'] == i][['x1','x2']].sum()/6, xytext=(-12,-4),textcoords='offset points',color='#404040',zorder=2)
        
    ax.scatter(dnodesnorm.x1,dnodesnorm.x2, s=20, c='#F0E68C',edgecolors='none', marker='o', zorder=1)   
    
    ax.set_xlim(-.05,1.05)
    ax.set_ylim(-.1,1.1)
    ax.set_frame_on(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.show()

    return()
  
colorRGB=[242,74,87, 124,218,198, 2,114,195]
plotallocation(7,2)
  
