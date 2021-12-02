# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 11:26:12 2020

@author: tajda
"""
'''
PLOTS THE FULL DISTRIBUTION OF NUMBER OF CASES AVERTED BY ONE PERSON TAKING ACTION TO AVOID INFECTION
'''

import numpy as np
import math as mt
import psutil as ps
import matplotlib.pyplot as plt

plt.close('all')

def Transmit_Set(R,psuc,errbnd,mu,sizemax, M):
    '''
    
    #---------------------------------Inputs----------------------------------#
    R: The reproduction number. Note that this is tied to mu
    psuc: The probability that a transmission will succeed (dependent on R)
    errbnd: The acceptable error.
    mu: The offspring distribution PGF. Note that the function must be able to
        calculate for a vector in one go.
    genrange: The generation range for which the calculation is needed. The 
        first generation is called generation 0.
    sizemax: The maximum transmission set size for which the probability is
        calulated. It should be considerably less than M (an order of 
        magnitude).
    
    #---------------------------Generated variables---------------------------#
    #gmax: The number of gens needed beyond the final gen considered
    #M: The number of points used in the estimate
    #m: a vector of numbers 0 to M-1 - used to generate z
    #r: the radius of the integral estimated in the approximation
    #z: the points used in the integral estimate
    #omega: this variable contains the Omega(z) for all z values. It changes 
        #depending on the generation of Omega.
    
    #---------------------Variables that can be adjusted----------------------#
    #print('beginning function')
    '''

    r=0.99 #This number can probably be anything less than 1.
    lenR=len(R)
    R=R.reshape(1,lenR,1)
    psuc=psuc.reshape(1,lenR,1)         
    if sizemax>=0.2*M:
        M=mt.ceil(sizemax*36/7)
        print('\n The number of points in the approximation has been increased due to the need to calculate large transmission set size(s) being considered. This may lead to a very long calculation time.')
    if M*lenR*(sizemax+1)*16>0.99*ps.virtual_memory().available:
        raise Exception('Calculation may exceed available memory')
    #Calculating the maximum generation
    gmax=mt.ceil(np.max(np.log(errbnd*(1-R*psuc))/np.log(R*psuc)))
    
    #Set-up for finding coefficients
    m=np.arange(M)
    z=(r*np.exp(2*mt.pi*1j*m/M)).reshape(1,1,M)
    omegaz=z
    for i in range(gmax):
        omegaz=z*mu(1-psuc+psuc*omegaz,R)
    size=np.arange(sizemax+1).reshape(sizemax+1,1,1)
    denoms=z**size
    terms=omegaz/denoms
    probs=1/M*np.sum(terms,axis=2).real
    return probs

def P_Saved_Poisson(R,errbnd):
    beta=np.zeros(len(R))
    beta[R<=1]=1
    while np.max(np.logical_and(np.exp(R[R>1]*(beta[R>1]+errbnd-1))-beta[R>1]-errbnd>0,beta[R>1]+errbnd<1))==True:
        beta[R>1]=np.exp(R[R>1]*(beta[R>1]-1))
    return beta

def mu_p(x,Rc):
    return np.exp(Rc*(x-1))

theta=0.16
def mu_nb(x,Rc):
    r=theta
    var=Rc+1/r*Rc**2
    p=(var-Rc)/var
    q=1-p
    return (q/(1-p*x))**r



Rc=np.array([0.75, 0.95, 1, 1.1, 1.5, 2.0, 2.5, 3.5])

decplaces=10 #Note that if this is more than ~8, you don't seem to get a higher level of accuracy anyway, because of rounding error.
errbnd=0.5*10**(-decplaces)
ps_p=P_Saved_Poisson(Rc,errbnd)
sizemax=1000
M=5003 #It's good to use a prime here, which should be more than 5 times the maximum size.

TSinfected_P=Transmit_Set(Rc,ps_p,errbnd,mu_p,sizemax,M)
TSinfected_NB=Transmit_Set(Rc,ps_p,errbnd,mu_nb,sizemax,M)

TStotal_P=(1-ps_p)*TSinfected_P
TStotal_NB=(1-ps_p)*TSinfected_NB

TStotal_P[0,:]=ps_p
TStotal_NB[0,:]=ps_p

'''
Plots
'''

def colour(val):
    c=np.zeros(3)
    if Rc[val]<=1:
        x=0.5*val/len(Rc[Rc<=1])
    else:
        x=0.5+0.5*(val-len(Rc[Rc<=1]))/len(Rc[Rc>1])
    if x<=1/3:
        c[1]=1-(1-3*x)**2
        c[2]=1-(3*x)**2
    elif x<=1/2:
        c[0]=1-(3-6*x)**2
        c[1]=1
    else:
        c[0]=1
        c[1]=2-2*x
    return [0.1,0,0]+[0.9,0.9,1]*c

numcases = np.arange(sizemax)

fig = plt.figure(figsize = (5.8, 7.5))
ax = fig.add_subplot(111)
ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

plt.gcf().subplots_adjust(top = 0.92, bottom = 0.09, left = 0.14, right = 0.75)
for i1 in range(len(Rc)):
    ax0.semilogy(numcases, TSinfected_P[1:,i1], '.-', label=r"$\mathcal{R}_c=" + str(Rc[i1]) + "$", c=colour(i1))
    ax1.semilogy(numcases, TSinfected_NB[1:,i1], '.-', label=r"$\mathcal{R}_c=" + str(Rc[i1]) + "$", c=colour(i1))

ax0.set_xscale('symlog')
ax1.set_xscale('symlog')

# ax.set_xlabel('Number of averted cases due to one \n individual isolating after infection')
#ax.set_xlabel('Number of cases averted by one infected individual isolating')
ax.set_xlabel('Number of infections averted')

ax.set_ylabel('Probability')
ax.yaxis.labelpad = 15

ax0.set_xlim([0, sizemax])
ax1.set_xlim([0, sizemax]) 
ax0.set_ylim([1e-8, 1])
ax1.set_ylim([1e-8, 1])

ax0.set_title('Poisson offspring distribution')
ax1.set_title('Negative binomial offspring distribution')
ax1.legend(bbox_to_anchor=(1.03, 1.1), loc='center left')

plt.savefig('distsymloglog.pdf')