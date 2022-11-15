# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:43:21 2020

@author: tajda
"""
'''
PLOTS 3 HEAT MAPS SHOWING DIFFERENT STATISTICS AGAINST STARTING REPRODUCTION NUMBER AND PROPORTION OF PEOPLE WHO ISOLATE:
1. NUMBER OF INFECTIONS PREVENTED PER ISOLATING INDIVIDUAL
2. NUMBER OF INFECTIONS PREVENTED BY ONE MORE PERSON ISOLATING
3. PROPORTION OF TOTAL POPULATION ISOLATING (THIS ONE ISN'T IN THE PAPER)  [no longer does this (commented out at bottom)]
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl

plt.close('all')

'''
This function calculates the probability that a random person is infected in an
epidemic, assuming everyone is equally susceptible.
This is the same as the proportion of people infected in the large population
limit
This is performed for a set of reproduction numbers, supplied as the vector R0.
'''
def P_Infected(R0, errbnd):
    beta = np.zeros_like(R0)
    beta[R0 <= 1] = 1
    R0g1 = R0[R0 > 1]
    betag1 = beta[R0 > 1]
    count = 0
    while np.max(np.logical_and(np.exp(R0g1*(betag1 + errbnd - 1)) - betag1 - errbnd > 0, betag1 + errbnd < 1)) == True:
        betag1 = (R0g1*betag1-1)*np.exp(R0g1*(betag1 - 1))/(R0g1*np.exp(R0g1*(betag1-1))-1)
        beta[R0 > 1] = betag1
        count += 1
    return 1-beta

'''
This function makes a heatmap from a set of data.
xvals/yvals: Vectors giving the x- and y-values (dependent variables)
data: A 2D array containing all dependent variable values
title: The figure title
xlab/ylab: x- and y-axis labels
cbarlabel: A label for the colour bar
'''
def graphmaker(ax,xvals, yvals, data, title, xlab, ylab, cbarlabel, data_min, data_max, levels,**kwargs):
    x,y=np.meshgrid(xvals,yvals)
    
    if 'norm' in kwargs and kwargs.get('norm')=='log':
        norm=cl.LogNorm(vmin=data_min, vmax=data_max)
    else:
        norm=cl.Normalize(vmin=data_min,vmax=data_max)
    im=ax.pcolormesh(x,y,data,cmap="Spectral_r",norm=norm)
    cbar=plt.colorbar(im, extend='max', ax=ax)
    
            
    contours = ax.contour(x, y, data, levels, colors='k', linewidths=0.75)
    plt.clabel(contours, rightside_up=True, fmt='%g')
    
    cbar.set_label(cbarlabel, FontSize=12)
    #plt.suptitle(title, fontsize=11)
    ax.ticklabel_format(style='sci',scilimits=(-2,2))
    ax.set_xlabel(xlab, FontSize=12)
    ax.set_ylabel(ylab, FontSize=12)

    ax.set_xlim([np.min(xvals),np.max(xvals)])
    ax.set_ylim([np.min(yvals),np.max(yvals)])
    #ax.subplots_adjust(left=0.1,right=0.95,bottom=0.15,wspace=0.25)

'''
Initial parameters
'''
res = 600
decplaces = 10 #The number of decimal placed of accuracy. 
             #Note that actual accuracy does not seem to go beyond 8 decimal 
             #places, due to rounding error.
errbnd = 0.5*10**(-decplaces) #The error bound.

'''
R0: The basic reproduction number (reproduction number with no herd immunity 
    and no interventions)
c:  The proportion of transmissions blocked by control measures
Rc: The reproduction number under control (reproduction number when control 
    measures block a proportion c of transmissions)
'''
R0 = np.linspace(1, 4, res) 
c = np.linspace(0, 1, res)
Rc = R0*(1-c.reshape(res,1))


fig = plt.figure(figsize = (14.5, 4))
axL = fig.add_subplot(121)
axR = fig.add_subplot(122)
plt.subplots_adjust(wspace=0.25)


prop_infected_c = P_Infected(Rc,errbnd)
prop_infected_none = prop_infected_c[0,:]
prop_saved = prop_infected_none - prop_infected_c

num_saved_per_person = prop_saved/(prop_infected_c*c)

graphmaker(axL, R0, c, num_saved_per_person, '', 'Starting reproduction number', 'Proportion of infected people who isolate', 'Number of infections averted \n per isolating individual', 1e-3, 1e3, [1e-2, 1e-1, 1, 10, 1e2], norm='log')

#plt.savefig('ManyPeopleIsolating.png')

num_saved_by_one_more_isolating = R0*(1 - prop_infected_c)/(1 - Rc*(1 - prop_infected_c))

graphmaker(axR, R0, c, num_saved_by_one_more_isolating, '', 'Starting reproduction number', 'Proportion of infected people who isolate', 'Number of infections averted by \n one more individual isolating', 1e-1, 1e3, [0.1, 0.3, 1, 3, 10, 100], norm='log')

plt.savefig('Isolating.png', bbox_inches='tight')

#prop_isolating = prop_infected_c*c.reshape(res, 1)

#graphmaker(R0, c, prop_isolating, '', 'Starting reproduction number', 'Proportion of infected people who isolate', 'Proportion of total population isolating', 0, 0.41, [0, 0.1, 0.2, 0.3, 0.4])

#plt.savefig('PropIsolating.png')
