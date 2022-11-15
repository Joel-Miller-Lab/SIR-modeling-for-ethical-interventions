# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 18:18:13 2020

@author: tajda
"""
'''
CALCULATES AND PLOTS THE TWO PLOTS RELATED TO "MARGINAL BENEFIT" OF ONE PERSON ISOLATING, 
DEPENDING ON THE ACTIONS OF OTHER INFECTED INDIVIDUALS.
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

Rc = np.array([1.1,1.5,2,2.5,3.5])
lenRc=len(Rc)
decplaces = 10
errbnd = 10**(-decplaces)
maxit = 1000
i = 0

theta = 0.16

pep = np.ones(lenRc).reshape(1,lenRc) #probability of epidemic for Poisson
while np.min(1 - np.exp(-Rc*(pep - errbnd)) - pep + errbnd) < 0 and np.max(pep - errbnd) > 0 and i < maxit:
    pep = ((Rc*pep + 1)*np.exp(-Rc*pep) - 1)/(Rc*np.exp(-Rc*pep) - 1)
    i += 1

r, q = convert_R_to_numpy_params(Rc, theta)
p = 1 - q
i = 0

if i == maxit:
    raise Exception('Calculation took too long. (Poisson)')

penb = np.ones(lenRc).reshape(1,lenRc) #Probability of Epidemic for Negative Binomial
while np.min(1 - (q/(q + p*(penb - errbnd)))**r - penb + errbnd) < 0 and np.max(penb - errbnd) > 0 and i < maxit:
    penb = penb - (1 - (q/(q + p*penb))**r - penb)/(r*q**r*p/(q + p*penb)**(r + 1) - 1)
    i += 1

if i == maxit:
    raise Exception('Calculation took too long. (Negative Binomial)')

maxE=50
E=(np.arange(maxE+1)).reshape(maxE+1,1)

Benefit_p = pep*(1 - pep)**(E)
Benefit_nb = penb*(1 - penb)**(E)

fig = plt.figure(figsize = (11, 5.8))
axL = fig.add_subplot(121)
axR = fig.add_subplot(122)

ax0 = fig.add_subplot(221)
ax1 = fig.add_subplot(223)
ax0b = fig.add_subplot(222)
ax1b = fig.add_subplot(224)
plt.subplots_adjust(wspace=1.0,hspace=0.3)


axL.spines['top'].set_color('none')
axL.spines['bottom'].set_color('none')
axL.spines['left'].set_color('none')
axL.spines['right'].set_color('none')
axL.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

plt.gcf().subplots_adjust(top = 0.92, bottom = 0.09, left = 0.14, right = 0.77)

for i1 in range(lenRc):
    ax0.semilogy(E, Benefit_p[:,i1])
    ax1.semilogy(E, Benefit_nb[:,i1], label = r'$\mathcal{R}_c=' + str(Rc[i1]) +'$')

axL.set_xlabel('Number of other (non-isolating) introductions', FontSize=11)
axL.set_ylabel('Reduction in epidemic probability from one infected individual isolating', FontSize=11)
axL.yaxis.labelpad = 15
ax0.set_xlim([0, maxE])
ax1.set_xlim([0, maxE])
ax0.set_ylim([1e-10,2])
ax1.set_ylim(bottom = 1e-4)

ax0.set_title('Poisson offspring distribution', FontSize=11)
ax1.set_title('Negative binomial offspring distribution', FontSize=11)
ax1.legend(bbox_to_anchor=(1.01, 1.58), loc='center left')

#plt.savefig('MarginalBenefit.png')

Tot_infected = 50
Isolating = np.arange(Tot_infected + 1).reshape(Tot_infected+1,1)

Av_reduction_p=1/Isolating[1:]*(1 - (1 - pep)**Isolating[1:])*(1 - pep)**(Tot_infected - Isolating[1:])
Next_margben_p=pep*(1 - pep)**(Tot_infected - Isolating[:Tot_infected] - 1)

Av_reduction_nb=1/Isolating[1:]*(1 - (1 - penb)**Isolating[1:])*(1 - penb)**(Tot_infected - Isolating[1:])
Next_margben_nb=penb*(1 - penb)**(Tot_infected - Isolating[:Tot_infected] - 1)

Tot_benefit_p = (1-pep)**(Tot_infected-Isolating) - (1-pep)**(Tot_infected)
Tot_benefit_nb = (1-penb)**(Tot_infected-Isolating)- (1-penb)**(Tot_infected)



axR.spines['top'].set_color('none')
axR.spines['bottom'].set_color('none')
axR.spines['left'].set_color('none')
axR.spines['right'].set_color('none')
axR.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

plt.gcf().subplots_adjust(top = 0.92, bottom = 0.09, left = 0.14, right = 0.77)
for i1 in range(lenRc):
    if i1 == 0:
        ax0b.semilogy(Isolating, Tot_benefit_p[:,i1], 'C' + str(i1),label='Total \nbenefit')
        ax0b.semilogy(Isolating[1:], Av_reduction_p[:,i1], 'C' + str(i1) + '--', label = 'Average \nbenefit')
        ax0b.semilogy(Isolating[:Tot_infected], Next_margben_p[:,i1], 'C' + str(i1)+ '-.', label = 'Marginal \nbenefit')
    else:
        ax0b.semilogy(Isolating[1:], Av_reduction_p[:,i1], 'C' + str(i1) + '--')
        ax0b.semilogy(Isolating[:Tot_infected], Next_margben_p[:,i1], 'C' + str(i1)+ '-.')
        ax0b.semilogy(Isolating, Tot_benefit_p[:,i1], 'C' + str(i1))
    ax1b.semilogy(Isolating[1:], Av_reduction_nb[:,i1], 'C' + str(i1) + '--')
    ax1b.semilogy(Isolating[:Tot_infected], Next_margben_nb[:,i1], 'C' + str(i1) + '-.')
    ax1b.semilogy(Isolating, Tot_benefit_nb[:,i1], 'C' + str(i1), label = r'$\mathcal{R}_c=' + str(Rc[i1]) +'$')

axR.set_xlabel('Number isolating (of 50 introductions)', FontSize=11)
axR.set_ylabel('Benefit', FontSize=11)
axR.yaxis.labelpad = 15

ax0b.set_ylim([1e-3,1])
ax1b.set_ylim([1e-3,1])

#ax0.set_xlim([Tot_infected-10, Tot_infected])
ax0b.set_xlim([35, Tot_infected])
ax1b.set_xlim([1, Tot_infected])
ax0b.set_title('Poisson offspring distribution', FontSize=11)
ax1b.set_title('Negative binomial offspring distribution', FontSize=11)
ax1b.legend(bbox_to_anchor=(1.01, 1.9), loc='upper left')
ax0b.legend(bbox_to_anchor=(1.01,-0.5), loc='upper left')

plt.savefig('marg' + str(Tot_infected) + '.png', bbox_inches='tight')
