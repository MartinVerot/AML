#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ionization energy without electron repulsion.

Informations
------------
Author : Martin Vérot  from the ENS de Lyon, France
Licence : Creative Commons CC-BY-NC-SA 4.0

"""

# Importation of libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import widgets
def Sintegral(zeta,R):
    return (1+zeta*R+zeta**2*R**2/3.)*np.exp(-zeta*R)

def Iintegral(zeta,R):
    return 1/R * (1-(1+zeta*R)*np.exp(-2*zeta*R))

def Jintegral(zeta,R):
    return zeta*(1+zeta*R)*np.exp(-zeta*R)

def Eg(zeta,R):
    return -zeta**2/2.+1/R+(zeta*(zeta-1)-Iintegral(zeta,R)+(zeta-2.)*Jintegral(zeta,R))/(1+Sintegral(zeta,R))

def Eg2(x):
    return Eg(x[0],x[1])

def Eu(zeta,R):
    return -zeta**2/2.+1/R+(zeta*(zeta-1)-Iintegral(zeta,R)-(zeta-2.)*Jintegral(zeta,R))/(1-Sintegral(zeta,R))

def createGrid(xmin,xmax,nx,ymin,ymax,ny):
    x, dx = np.linspace(xmin, xmax, nx, retstep = True)
    y, dy = np.linspace(ymin, ymax, ny, retstep = True)
    return np.meshgrid(x, y, indexing='xy'), (dx,dy)

# Main program
if __name__ == "__main__":

    #figure
    fig = plt.figure(figsize=(8,8))
    gs = fig.add_gridspec(2, 1,  left=0.08, right=0.95, bottom=0.05, top=0.95, wspace=0.18, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],projection='3d', proj_type='ortho')
    #ax3 = fig.add_subplot(gs[1,0])
    #ax4 = fig.add_subplot(gs[1,1])
    xs = np.linspace(0.0001,4,500)
    ax1.set_xlim(0.,4.)
    ax1.set_ylim(-0.7,0.5)
    grid, delta = createGrid(0, 2, 250,0.0001, 4, 250) 
    zetas,Rs = grid    
    Z = Eg(zetas,Rs)
    Z[Z>0]=0
    cutR=Rs[0::,0]
    #ax2.plot(,Eg(),label='$E_g$')
    #ax2.plot(xs,Eu(xs,0),label='$E_u$')
    ax2.set_xlabel('$\lambda$')
    ax2.set_ylabel('$R$')
    ax2.set_zlabel('$E$')
    #ax2.set_ylabel('$E$')
    ax2.set_zlim(np.min(Eg(zetas,Rs))-0.01,0.1)
    ax2.set_xlim(0.,2.)
    ax2.set_ylim(0.,4.)
   
    #Minimum global
    minglobal = minimize(Eg2,[1.,2.])

    #tracé du minimum le long de l'abscisse curviligne
    lines={}
    lines['surface'] = ax2.plot_surface(zetas, Rs, Z,ccount=100,rcount=100, cmap=cm.coolwarm, vmin=np.min(Eg(zetas,Rs)), vmax=0,zorder=-1,alpha=.5)#, linewidth=0, antialiased=False)
    #Tracé des énergies pour une valeur donnée de lambda
    lines2={}
    lines2['Egmin'],=ax1.plot(cutR,Eg(minglobal.x[0],cutR),label='optimum',linestyle='--', alpha=0.5)
    lines2['minfull'],=ax1.plot(minglobal.x[1],Eg(minglobal.x[0],minglobal.x[1]),label='R: {:.2f}, $\lambda$: {:.2f} E:{:.2f}'.format(minglobal.x[1],minglobal.x[0],minglobal.fun),linestyle='',marker='o', alpha=0.5)


 
    ax1.legend(loc='upper right')
    plt.show()
