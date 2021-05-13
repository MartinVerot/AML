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
import widgets

  
def BrancheP(K,omega_0,B,a):
    return 1/(constants.h*constants.c)*(constants.h*constants.c*omega_0+2*constants.h*constants.c*B*(K+1)-constants.h*constants.c*a*(K**2+4*K+3))

def BrancheR(K,omega_0,B,a):
    return 1/(constants.h*constants.c)*(constants.h*constants.c*omega_0-2*constants.h*constants.c*B*K-constants.h*constants.c*a*K*(K-2))
def pop(K,omega_0,B,a,T):
    return (2*K+1)*np.exp(-1/(constants.k*T)*(100*constants.h*constants.c*B*K*(K+1)-100*constants.h*constants.c*a*0.5*K*(K+1)))

def plot_data(omega_0,B,a,T):
    maxi = np.max([pop(xs,omega_0,B,a,T),pop(xs+1,omega_0,B,a,T)])
    lines['P'].set_data(BrancheP(xs,omega_0,B,a),pop(xs,omega_0,B,a,T))
    lines['R'].set_data(BrancheR(xs+1,omega_0,B,a),pop(xs+1, omega_0,B,a,T))
    lines['ref'].set_data([omega_0,omega_0],[0,0.25*maxi])
    ax1.set_xlim(2886.+300,2886.-300)
    ax1.set_ylim(0,maxi*1.1)
    for x in xs:
        textsR[int(x)].set_position((BrancheR(x+1,omega_0,B,a)-4, pop(x+1,omega_0,B,a,T) ))
        textsP[int(x)].set_position((BrancheP(x,omega_0,B,a)+10, pop(x,omega_0,B,a,T) ))
 
    
    
# Main program
if __name__ == "__main__":
    parameters = {
        'omega_0' : widgets.FloatSlider(value=2886., description='$\omega_0$', min=2700, max=3000),
        'B' : widgets.FloatSlider(value=10.75, description='$B$', min=0, max=20),
        'a' : widgets.FloatSlider(value=0., description='$a$', min=0, max=0.25),
        'T' : widgets.FloatSlider(value=300., description='$T$', min=100, max=500),        
    }

    fig = plt.figure(figsize=(8,8))
    gs = fig.add_gridspec(1, 1,  left=0.08, right=0.95, bottom=0.05, top=0.95, wspace=0.18, hspace=0.3)
    ax1 = fig.add_subplot(gs[0,0])
    #ax3 = fig.add_subplot(gs[1,0])
    #ax4 = fig.add_subplot(gs[1,1])
    xs = np.linspace(0,50,51)    
    
    lines={}    
    lines['P'], = ax1.plot([],[], marker='+', linestyle='',label = 'Branche P')
    lines['R'], = ax1.plot([],[], marker='+', linestyle='',label = 'Branche R')    
    lines['ref'], = ax1.plot([],[], marker='', linestyle='-',label = '$\omega_0$') 

    textsR=[]
    textsP=[]
    for x in xs:
        textsR.append(ax1.text(0, 0,  str('{}'.format(int(x)))))
        textsP.append(ax1.text(0, 0,  str('{}'.format(int(x+1)))))
    ax1.set_xlabel('E')
    ax1.legend(loc='upper right')
    param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.07, 0.96, 0.8, 0.03])
    #fig.suptitle(r"Énergie d'ionisation sans couplage e-e", fontsize=16)
    plt.show()
