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
   
def phi(r,R):
    return 1/np.sqrt(np.pi)*np.exp(-np.abs(r-R)) 

def Sintegral(R):
    return (1+R+R**2/3.)*np.exp(-R)

def Iintegral(R):
    return 1/R * (1-(1+R)*np.exp(-2*R))

def Jintegral(R):
    return (1+R)*np.exp(-R)

def Eg(R,E_0):
    return E_0+1/R-(Iintegral(R)+Jintegral(R))/(1+Sintegral(R))

def Eu(R,E_0):
    return E_0+1/R-(Iintegral(R)-Jintegral(R))/(1-Sintegral(R))

def Eg2(E_0,R):
    return E_0-(Iintegral(R)+Jintegral(R))/(1+Sintegral(R))

def Eu2(E_0,R):
    return E_0-(Iintegral(R)-Jintegral(R))/(1-Sintegral(R))


def Phi_g(r,R):
    return (1/np.sqrt(1+overlap(R))*(phi(r,0)+phi(r,R)))

def Phi_u(r,R):
    return (1/np.sqrt(1-overlap(R))*(phi(r,0)-phi(r,R)))
# Main program
if __name__ == "__main__":
    #figure
    fig = plt.figure(figsize=(8,8))
    gs = fig.add_gridspec(2, 1,  left=0.08, right=0.95, bottom=0.05, top=0.95, wspace=0.18, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    #ax3 = fig.add_subplot(gs[1,0])
    #ax4 = fig.add_subplot(gs[1,1])
    xs = np.linspace(0.0001,15,500)
    ax1.plot(xs,Sintegral(xs), label = 'S')
    ax1.plot(xs,Iintegral(xs), label = 'I')
    ax1.plot(xs,Jintegral(xs), label = 'J')
    ax1.set_xlabel('R')
    ax1.set_ylim(0,1.1)
    ax1.set_xlim(0,15)
    ax1.legend(loc='upper right')
    
    
    xs = np.linspace(0.5,10,500)
    ax2.plot(xs,Eg(xs,0),label='$E_g$')
    ax2.plot(xs,Eu(xs,0),label='$E_u$')
    ax2.set_xlabel('R')
    ax2.set_ylabel('$E$')
    ax2.set_ylim(np.min([Eg(xs,0),Eu(xs,0)])-0.01,0.1)
    #On trouve le minimum de E_g
    Emin = minimize(Eg,2.5,args=(0))    
    ax2.plot(Emin.x,[Emin.fun],label='R : {:.2f} E : {:.4f}'.format(Emin.x[0],Emin.fun),marker='o')
    #print(Emin)
    #print(Eg(Emin.x,0))
    ax2.legend(loc='upper right')

    #ax3.plot(xs,Eg2(0,xs),label='$E_g\'$')
    #ax3.plot(xs,Eu2(0,xs),label='$E_u\'$')
    #ax3.set_xlabel('R')
    #ax3.set_ylabel('$\Phi_\mathrm{u}$')
    #ax3.legend(loc='upper right')
    #ax3.set_ylim(np.min([Eg2(0,xs),Eu2(0,xs)]),0.4)
    #ax4.plot(xs,Phi_u(xs,3)**2)
    #ax4.set_xlabel('x')
    #ax4.set_ylabel('$|\Phi_\mathrm{u}|^2$')
    #ax4.set_ylim(0,0.5)
 
    #fig.suptitle(r"Énergie d'ionisation sans couplage e-e", fontsize=16)
    plt.show()
