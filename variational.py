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


def plot_data(lambdaa):
    energies = Eg(lambdaa,cutR)
    lines2['Eg'].set_data(cutR,energies)
    lines2['Eu'].set_data(cutR,Eu(lambdaa,cutR))
    lines2['min'].set_data(cutR[np.argmin(energies)],np.min(energies))

    energies2= Eg(lambdaa,cutR)
    energies2[energies2>0]=0.
    lines['plane'].set_data(lambdaa*np.ones_like(cutR),cutR)
    lines['plane'].set_3d_properties(energies2)     

# Main program
if __name__ == "__main__":
    parameters = {
        'lambdaa' : widgets.FloatSlider(value=1., description='$\lambda$', min=0.0001, max=2),
    }
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

    #ax2.plot(,Eg(),label='$E_g$')
    #ax2.plot(xs,Eu(xs,0),label='$E_u$')
    ax2.set_xlabel('$\lambda$')
    ax2.set_ylabel('$R$')
    ax2.set_zlabel('$E$')
    #ax2.set_ylabel('$E$')
    ax2.set_zlim(np.min(Eg(zetas,Rs))-0.01,0.1)
    ax2.set_xlim(0.,2.)
    ax2.set_ylim(0.,4.)
    #On trouve le minimum de E_g
    #Emin = minimize(Eg,2.5,args=(0))    
    zetaMins=[]
    EMins=[]
    cutR=Rs[0::,0]
    for R in cutR:
        Emin=minimize(Eg, x0 = 1. ,args=R)
        zetaMins.append(Emin.x[0])
        EMins.append(Emin.fun)
    EMins=np.asarray(EMins)
    EMins[EMins>0]=0

    #Minimum global
    minglobal = minimize(Eg2,[1.,2.])
    #print(minglobal)

    #tracé du minimum le long de l'abscisse curviligne
    lines={}
    lines['optimum'], = ax2.plot(zetaMins,Rs[0::,0],EMins, zorder=0,color='black')
    lines['surface'] = ax2.plot_surface(zetas, Rs, Z,ccount=100,rcount=100, cmap=cm.coolwarm, vmin=np.min(Eg(zetas,Rs)), vmax=0,zorder=-1,alpha=.5)#, linewidth=0, antialiased=False)
    lines['plane'], = ax2.plot([],[],[], alpha=.75,color='C0')#, linewidth=0, antialiased=False)
    #Tracé des énergies pour une valeur donnée de lambda
    lines2={}
    lines2['Eg'],=ax1.plot([],[],label='$E_\mathrm{g}$',color='C0')
    lines2['Eu'],=ax1.plot([],[],label='$E_\mathrm{u}$',color='C1')
    lines2['min'],=ax1.plot([],[], marker='o')
    lines2['Egmin'],=ax1.plot(cutR,Eg(minglobal.x[0],cutR),label='optimum',linestyle='--', alpha=0.5)
    lines2['minfull'],=ax1.plot(minglobal.x[1],Eg(minglobal.x[0],minglobal.x[1]),label='R: {:.2f}, $\lambda$: {:.2f} E:{:.2f}'.format(minglobal.x[1],minglobal.x[0],minglobal.fun),linestyle='',marker='o', alpha=0.5)
    lines2['Egmin2'],=ax1.plot(cutR,Eg(np.asarray(zetaMins),cutR),label='optimum2',color='black',linestyle='--', alpha=0.5)

    #box pour le choix des courbes
    choose_widget = widgets.make_choose_plot(lines, box=[0.005, 0.1, 0.12, 0.08])

    #box pour les paramètre
    param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.96, 0.4, 0.02])
    
    ax1.legend(loc='upper right')
    plt.show()
