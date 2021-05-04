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
   
def phi(r,R):
    return 1/np.sqrt(np.pi)*np.exp(-np.abs(r-R)) 

def overlap(R):
    return (1+R+R**2/3.)*np.exp(-R)

def Phi_g(r,R):
    return (1/np.sqrt(1+overlap(R))*(phi(r,0)+phi(r,R)))

def Phi_u(r,R):
    return (1/np.sqrt(1-overlap(R))*(phi(r,0)-phi(r,R)))
# Main program
if __name__ == "__main__":
    #figure
    fig = plt.figure(figsize=(10,6))
    gs = fig.add_gridspec(2, 2,  left=0.08, right=0.95, bottom=0.05, top=0.85, wspace=0.18, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    xs = np.linspace(-2,5,500) 
    ax1.plot(xs,Phi_g(xs,3))
    ax1.set_xlabel('x')
    ax1.set_ylabel('$\Phi_\mathrm{g}$')
    ax1.set_ylim(-0.7,0.7)
    
    ax2.plot(xs,Phi_u(xs,3))
    ax2.set_xlabel('x')
    ax2.set_ylabel('$\Phi_\mathrm{u}$')
    ax2.set_ylim(-0.7,0.7)

    ax3.plot(xs,Phi_g(xs,3)**2)
    ax3.set_xlabel('x')
    ax3.set_ylabel('$|\Phi_\mathrm{g}|^2$')
    ax3.set_ylim(0,0.5)
    
    ax4.plot(xs,Phi_u(xs,3)**2)
    ax4.set_xlabel('x')
    ax4.set_ylabel('$|\Phi_\mathrm{u}|^2$')
    ax4.set_ylim(0,0.5)
 
    #fig.suptitle(r"Énergie d'ionisation sans couplage e-e", fontsize=16)
    plt.show()
