#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ionization energy without electron repulsion.

Informations
------------
Author : Martin Vérot  from the ENS de Lyon, France, with the help of some scripts to create the buttons and sliders taken from https://github.com/araoux/python_agregation (written by P Cladé, A Raoux and F Levrier)
Licence : Creative Commons CC-BY-NC-SA 4.0

WARNING this program requires the widgets.py file to work
"""

# Importation of libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.integrate
import widgets
import matplotlib.patches as patches
#plt.rc('text', usetex=True)
#plt.rc('text.latex', preamble=r'\usepackage{stackrel}')




def plot_data(Z):
    #Maximum occupancy for each value of n
    Threshold=[0,2,8,18,32,50,72]
    CumSum = np.cumsum(Threshold)
    k=0
    EnIo = 0.
    Elements = ['', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    for n in range(1,len(Threshold)):
        #Finding the current ionization energy for the given value of Z
        if Z >CumSum[n-1]:
            EnIo = Z**2/n**2
    #changing the fill color of the rectangle
    for rect in rects:
        if rects.index(rect)+1==Z:
            rect.set_facecolor('red')
        elif rects.index(rect)+1<Z:
            rect.set_facecolor('#dddddd')
        else:
            rect.set_facecolor('white')
    #print(EnIo)
    #Updating the dots
    lines['Z'].set_data([Z],[EnIo])  
    lines['Z2'].set_data([Z],[-EnIo])
    #Updatind the vertical line
    lines['Z1'].set_data([Z,Z],[0,-30000])
    t1.set_position((Z+1, EnIo-30 ))
    t1.set_text(Elements[int(Z)])
    fig.canvas.draw_idle()
    
    
# Main program
if __name__ == "__main__":


    parameters = {
        'Z' : widgets.IntSlider(value=1, description='$Z$', min=1, max=118),
    }
    #Maximum occupancy for each value of n
    Threshold=[0,2,8,18,32,50,72]#,98,128]
    #Maximum number of electrons for n fixed.
    CumSum = np.cumsum(Threshold)
    
    
    fig = plt.figure(figsize=(10,6))
    #fig,axes=plt.subplots(2,2)
    gs = fig.add_gridspec(2, 2,  left=0.08, right=0.95, bottom=0.05, top=0.85, wspace=0.18, hspace=0.3)
    
    #Énergie des orbitales
    ax1 = fig.add_subplot(gs[0,0])#plt.subplot(2,2,1)
    ax1.set_ylim(-500,0)
    #Remplissage des couches de n fixé
    ax2 = fig.add_subplot(gs[0::,1])#plt.subplot(2,2,2)
    ax2.set_xlim(0,np.max(Threshold))
    ax2.set_ylim(0.5,len(Threshold)-0.5)
    ax2.set_title('Remplissage des couches')

    Zs = np.linspace(1,118,118)
    #print(Zs)
    #print(CumSum)
    
    #On trace les énergies orbitalaires
    for i in range(len(Threshold)):
        cutZ=Zs[Zs>CumSum[i]]
        #print('{}'.format(i))
        #print(cutZ)
        ax1.plot(cutZ,-cutZ**2/(i+1)**2,label='n = {:d}'.format(i+1),color='C{:d}'.format(i))
    
    #Énergie d'ionisation
    ax3 = fig.add_subplot(gs[1,0])#plt.subplot(2,2,3)
    Ei = np.zeros_like(Zs)
    for m in range(len(Zs)):
        for n in range(1,len(Threshold)):
            if Zs[m]>CumSum[n-1]:
                Ei[m]=Zs[m]**2/n**2
                #print('Z : {} n: {} Ei: {}'.format(Zs[m],n,Ei[m]))
    ax3.plot(Zs,Ei)
    #print(Ei)
    rects=[]
    #Tracé des rectangles
    for n in range(1,len(Threshold)):
        for j in range(Threshold[n]):
            if j < 2:
                edgecolor = '#4aae4a'
            elif j < 8:
                edgecolor = '#bcaed4'
            elif j < 18:
                edgecolor = '#fdc086'
            elif j < 32:
                edgecolor = '#295083'
            elif j<50:
                edgecolor = '#ffed6f'
            elif j<72:    
                edgecolor = 'black'
            #Plotting the filled rectangles for the current value of Z. The last electron is in red, the filled spinorbitals are in black and the empty ones are in white
            rects.append(ax2.add_patch(patches.Rectangle((j, n-0.5), 1,1, edgecolor = edgecolor, facecolor = 'white', fill=True)))

    ax1.set_xlim(min(Zs),max(Zs))
    ax3.set_xlim(min(Zs),max(Zs))
    ax1.set_xlabel('Z')
    ax1.set_ylabel('$E(n)$ (Ry)')
    ax3.set_xlabel('Z')
    ax3.set_ylabel('$E_I$ (Ry)')

    lines = {}
    #Dot for the ionization energy
    lines['Z'], = ax3.plot([],[],label='',color='C0',marker='o')
    #Dot for the orbital energy
    lines['Z2'], = ax1.plot([],[],label='',color='C0',marker='o')
    #Vertical line to show where is the Z value
    lines['Z1'], = ax1.plot([],[],label='',color='C0')
    t1 = ax3.text(0, 0,  str(''))
    ax1.legend(loc='upper right')

    fig.suptitle(r"Énergie d'ionisation sans couplage e-e", fontsize=16)
    param_widgets = widgets.make_param_widgets(parameters, plot_data, slider_box=[0.35, 0.91, 0.4, 0.02])


    plt.show()
