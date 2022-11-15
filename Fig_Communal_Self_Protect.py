# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:14:17 2020

@author: tajda & Joel Miller
"""

'''
CREATES 2 PLOTS
- FIRST
- SECOND
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl

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

def unvax_inf_prob(pre_v_R, v_f, its = 1000):
    #calculates the proportion of the entire population infected
    new_R = pre_v_R*(1-v_f)
    final_unvax_uninf_prop = np.zeros(new_R.shape)
    final_unvax_uninf_prop[new_R<1] = 1
    for ctr in range(its):
        #print(ctr)
        final_unvax_uninf_prop = np.exp(-new_R*(1-final_unvax_uninf_prop))
    return 1-final_unvax_uninf_prop

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
r_init = np.linspace(1, 4, 101)
vaxfrac=np.linspace(0, 1, 101)
original_fs= unvax_inf_prob(r_init, 0)

Orig_fs, VaxFrac = np.meshgrid(original_fs, vaxfrac) #there's probably a better way than
                                                     #meshgrid since I only want the final size
R_Init, VaxFrac = np.meshgrid(r_init, vaxfrac )
Unvax_final_size = unvax_inf_prob(R_Init, VaxFrac)    #this calculates the final proportion
                                        # in the population restricted to unvax

New_fs = Unvax_final_size*(1-VaxFrac)
impact = Orig_fs - New_fs
average_impact = impact/(VaxFrac+10**(-10))

R_final = R_Init * (1- VaxFrac - New_fs)

marginal_direct_vaccine_protection = Unvax_final_size #this is how much protection 1 more vaccine provides the
                                           #vaccinee
marginal_indirect_vaccine_protection = marginal_direct_vaccine_protection*(R_final/(1-R_final))
            #This is how much indirect protection comes from vaccinee [after accounting for P(infection)]
marginal_total_vaccine_protection = marginal_indirect_vaccine_protection+marginal_direct_vaccine_protection
'''
PLOTTING GRAPHS
'''


fig = plt.figure(figsize = (14.5, 4))
axL = fig.add_subplot(121)
axR = fig.add_subplot(122)
plt.subplots_adjust(wspace=0.25)


#fig, ax = plt.subplots(figsize = (5.8, 3.75))
data_min = 0#average_impact[average_impact>0].min()
data_max = 2#average_impact.max()
norm=cl.Normalize(vmin=data_min,vmax=data_max)
im = axR.pcolormesh(R_Init, VaxFrac, average_impact, cmap="Spectral_r", norm=norm)
cbar = fig.colorbar(im, extend='max', ax=axR)
contours = axR.contour(R_Init, VaxFrac, average_impact, np.linspace(0,2,11), colors='k', linewidths=0.75)
plt.clabel(contours, rightside_up=True, manual=False, fmt='%g')

#cbar.set_label(cbarlabel)
cbar.ax.set_ylabel('Average number of infections averted\n per individual avoiding infection', FontSize=12)#cbar.set_label('B')
#ax.ticklabel_format(style='sci', scilimits=(-2, 2))

axR.set_xlim([r_init.min(), r_init.max()])
axR.set_ylim([vaxfrac.min(), vaxfrac.max()])
#axR.set_subplots_adjust(left=0.1, right=0.95, bottom=0.15, wspace=0.25)
axR.set_xlabel('Starting Reproduction Number', FontSize=12)
axR.set_ylabel('Fraction avoiding infection', FontSize=12)


data_min = marginal_total_vaccine_protection.min()
data_max = marginal_total_vaccine_protection.max()
norm=cl.Normalize(vmin=data_min,vmax=data_max)
im = axL.pcolormesh(R_Init, VaxFrac, marginal_total_vaccine_protection, cmap="Spectral_r", norm=norm)
cbar = fig.colorbar(im, extend='max', ax=axL)
contours = axL.contour(R_Init, VaxFrac, marginal_total_vaccine_protection, np.linspace(0,2,11), colors='k', linewidths=0.75)
plt.clabel(contours, rightside_up=True, manual=False, fmt='%g')
cbar.ax.set_ylabel('Expected number of infections averted\n by one more individual avoiding infection', FontSize=12)#cbar.set_label('B')

#cbar.set_label(cbarlabel)
#fig.suptitle(title, fontsize=12)
#plt.ticklabel_format(style='sci', scilimits=(-2, 2))
#cbar.set_label("A")
#plt.xlabel(xlab, FontSize=12)
#plt.ylabel(ylab, FontSize=12)

axL.set_xlim([r_init.min(), r_init.max()])
axL.set_ylim([vaxfrac.min(), vaxfrac.max()])
#plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, wspace=0.25)
axL.set_xlabel('Starting Reproduction Number', FontSize=12)
axL.set_ylabel('Fraction avoiding infection', FontSize=12)
plt.savefig('tmp2.png', bbox_inches='tight')

# plt.gcf().subplots_adjust(top = 0.92, bottom = 0.15, left = 0.15, right = 0.95)
# ax.plot(R0, E_infected)
# ax.set_xlabel('Starting reproduction number')
# ax.set_ylabel('Fraction taking protective measures')
# ax.set_xlim([np.min(R0), np.max(R0)])
# ax.set_ylim([0.05,100])
# ax.set_yticks([1,10,100])
# ax.set_yticklabels(['1', '10', '100'])
# plt.savefig('tmp1.png')
#
# fig, ax = plt.subplots(figsize = (5.8, 3.75))
# plt.gcf().subplots_adjust(top = 0.92, bottom = 0.15, left = 0.1, right = 0.95)
# ax.plot(R0[R0 > 1], E_total[R0 > 1], label = 'Expected number saved by avoiding infection')
# ax.plot(R0[R0 > 1], 1 - ps[R0 > 1], label = 'Probability of being infected if not avoiding infection', FontSize=12)
# ax.set_xlabel('Reproduction number')
# ax.set_xlim([np.min(R0[R0 > 1]), np.max(R0)])
# ax.set_ylim([0, 2])
# ax.legend(loc = 'lower right')
# plt.savefig('ExpNumSaved.png')


#plt.show()


