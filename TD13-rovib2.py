#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Integrals for H_2^+ S, I, J, E

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
import scipy.constants as constants
  
def BrancheP(K,omega_0,B,a):
    return 1/(constants.h*constants.c)*(constants.h*constants.c*omega_0+2*constants.h*constants.c*B*(K+1)-constants.h*constants.c*a*(K**2+4*K+3))

def BrancheR(K,omega_0,B,a):
    return 1/(constants.h*constants.c)*(constants.h*constants.c*omega_0-2*constants.h*constants.c*B*K-constants.h*constants.c*a*K*(K-2))

# Main program
if __name__ == "__main__":
    #figure
    omega_0 =2886.
    B = 10.75
    a = 0.5
    fig = plt.figure(figsize=(8,8))
    gs = fig.add_gridspec(1, 1,  left=0.08, right=0.95, bottom=0.05, top=0.95, wspace=0.18, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    #ax3 = fig.add_subplot(gs[1,0])
    #ax4 = fig.add_subplot(gs[1,1])
    xs = np.linspace(0,12,13)
    ax1.plot(BrancheP(xs,omega_0,B,a),xs, marker='+', linestyle='',label = '$f(K)$')
    #xs = np.linspace(1,13,13)
    ax1.plot(BrancheR(xs,omega_0,B,a),xs, marker='+', linestyle='',label = '$f(K)$')
    ax1.set_xlabel('E')
    ax1.legend(loc='upper right')

    #fig.suptitle(r"Énergie d'ionisation sans couplage e-e", fontsize=16)
    plt.show()
