# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:14:17 2020

@author: tajda
"""

'''
CREATES 2 PLOTS
- FIRST SHOWS THE EXPECTED NUMBER OF CASES AVERTED BY AN INFECTED PERSON ISOLATING EFFECTIVELY
- SECOND SHOWS EXPECTED NUMBER OF CASES AVERTED BY SOMEONE TAKING SUFFICIENT ACTION TO PREVENT INFECTION
'''

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

'''
This function calculates the probability that a random person avoids infection
in an epidemic, assuming everyone is equally susceptible.
If a certain person is infected and passes the disease to another person, this
is also the probability that the second person would have been saved (would not
have been infected at all) if the first person isolated. 
This is performed for a set of reproduction numbers, supplied as the vector R0.
'''
def P_Saved(R0, errbnd):
    beta = np.zeros(len(R0))
    beta[R0 <= 1] = 1
    R0g1 = R0[R0 > 1]
    betag1 = beta[R0 > 1]
    while np.max(np.logical_and(np.exp(R0g1*(betag1 + errbnd - 1)) - betag1 - errbnd > 0, betag1 + errbnd < 1)) == True:
        betag1 = np.exp(R0g1*(betag1 - 1))
    beta[R0 > 1] = betag1
    return beta

'''
Initial parameters
'''
res = 1000 
R0 = (np.arange(res) + 1)/res*4 #A vector containing the reproduction numbers

decplaces = 10 #The number of decimal placed of accuracy. 
             #Note that actual accuracy does not seem to go beyond 8 decimal 
             #places, due to rounding error.
errbnd = 0.5*10**(-decplaces) #The error bound.

'''
CALCULATIONS
We calculate the number of further cases an infected individual will be
solely responsible for, i.e. the number of people who would not have become
infected if he had self-isolated and avoided infection.
We then calculate the "total risk" by multiplying this by the probability of
becoming infected. This gives the expectation value for the number of cases
that any random non-isolating individual is responsible for, including both 
people who become infected and people who avoid it.
Note that the individual is included in these expectation values.
'''
ps = P_Saved(R0, errbnd) #The probability of avoiding infection
E_infected = R0*ps/(1 - R0*ps) #Expected number caused by an infected person
E_total=(1 - ps)/(1 - R0*ps) #Total risk

'''
PLOTTING GRAPHS
'''

fig, ax = plt.subplots(figsize = (5.8, 3.75))
plt.gcf().subplots_adjust(top = 0.92, bottom = 0.15, left = 0.15, right = 0.95)
ax.semilogy(R0, E_infected)
ax.set_xlabel('Reproduction number')
ax.set_ylabel('Expected number of averted infections')
ax.set_xlim([np.min(R0), np.max(R0)])
ax.set_ylim([0.05,100])
ax.set_yticks([1,10,100])
ax.set_yticklabels(['1', '10', '100'])
plt.savefig('ExpNumSavedifinfected.png')

fig, ax = plt.subplots(figsize = (5.8, 3.75))
plt.gcf().subplots_adjust(top = 0.92, bottom = 0.15, left = 0.1, right = 0.95)
ax.plot(R0[R0 > 1], E_total[R0 > 1], label = 'Expected total infections averted')
ax.plot(R0[R0 > 1], 1 - ps[R0 > 1], label = 'Probability of being infected if not avoiding infection')
ax.fill_between(R0[R0 > 1], E_total[R0 > 1], 1 - ps[R0 > 1], color = 'yellow')
ax.set_ylabel('Expected Benefit')
ax.set_xlabel('Reproduction number')
ax.set_xlim([np.min(R0[R0 > 1]), np.max(R0)])
ax.set_ylim([0, 2])
ax.legend(loc = 'lower right')
plt.savefig('ExpNumSaved.png')
