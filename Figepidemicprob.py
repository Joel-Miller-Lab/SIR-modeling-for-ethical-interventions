# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 13:39:28 2020

@author: tajda
"""

'''
CALCULATES AND PLOTS THE PROBABILITY THAT ONE IMPORTED CASE CAUSES AN EPIDEMIC.
'''


import matplotlib.pyplot as plt
import numpy as np

def convert_R_to_numpy_params(mu, theta):
    """
    Convert mean/dispersion parameterization of a negative binomial to the ones numpy supports

    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
    From https://stackoverflow.com/a/47406400/2966723
    """
    r = theta
    var = mu + 1 / r * mu ** 2
    p = (var - mu) / var
    return r, 1 - p

plt.close('all')

res = 1000
R0 = (np.arange(res) + 1)/res*3 + 1
decplaces = 10
errbnd = 10**(-decplaces)
maxit = 1000
i = 0

theta = 0.16

pep = np.ones(res)
while np.min(1 - np.exp(-R0*(pep - errbnd)) - pep + errbnd) < 0 and np.max(pep - errbnd) > 0 and i < maxit:
    pep = ((R0*pep + 1)*np.exp(-R0*pep) - 1)/(R0*np.exp(-R0*pep) - 1)
    i += 1

r, q = convert_R_to_numpy_params(R0, theta)
p = 1 - q
i = 0

if i == maxit:
    raise Exception('Calculation took too long. (Poisson)')

penb = np.ones(res)
while np.min(1 - (q/(q + p*(penb - errbnd)))**r - penb + errbnd) < 0 and np.max(penb - errbnd) > 0 and i < maxit:
    penb = penb - (1 - (q/(q + p*penb))**r - penb)/(r*q**r*p/(q + p*penb)**(r + 1) - 1)
    i += 1

if i == maxit:
    raise Exception('Calculation took too long. (Negative Binomial)')


plt.figure(figsize = (5.8,3.75))
plt.plot(R0, pep, label = 'Poisson')
plt.plot(R0, penb, label = 'Negative Binomial')
plt.gcf().subplots_adjust(bottom = 0.15)
plt.gcf().subplots_adjust(left = 0.15)
plt.legend()
plt.xlabel('Reproduction number')
plt.ylabel('Probability of causing an epidemic')
name = 'Probofepidemic.png'
plt.xlim([1, np.max(R0)])
plt.ylim([0, 1])
plt.savefig(name)