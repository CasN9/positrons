# -*- coding: utf-8 -*-
# Organisation is always a good thing. Make each step a function?
# May not be the best thing to do in python... (Matlab tip). 

# Caspian Nicholls, u1027945
# Australian National University

# Reference for vibrational frequencies from http://cccbdb.nist.gov/alldata1.asp
import math
import pylab
import matplotlib.pyplot as plt
#import numpy.random as rnd
from scipy.integrate import trapz, simps
import csv

## Graphs to use: 
## Graph of binding energy against dipole polarizability
## Graph of amplitude of C-H stretch peak against number of carbons
## Graphs of number against size or mass
## Graph of Z_eff against positron energy
#
## Given a positron energy and PAH size, obtain dipole polarizability and number
## of carbons, then binding energy and C-H stretch peak amplitude. 
## Goal: calculate positron lifetime for different alpha values.
#
#### Constants ###
## The first three are not realistic values. 
#density = 1  # number density
#sigma = 1    # cross-section
#vel = 1      # velocity
#rate = sigma*density*vel  # collision rate
#alpha = 3.0/2.0  # dipole polarizability (3/2, 7/2, 19/12, 5/2 etc.)
#char_width = 10e-6  # characteristic width of resonance, Gamma_v
#exp_width = 10e-2  # experimental width of resonance
#
#
def energy_dist(Z_eff, gamma, E, E_0):
    '''Gaussian with variable x - E_0 is incident energy, E is positron energy.'''
    coeff = Z_eff / (gamma * math.sqrt(2 * math.pi))
    power1 = - (E - E_0)**2
    power2 = 2* (gamma)**2
    f_E = coeff * math.exp(power1/ power2)
    return f_E

#
## Plotting the size distribution (radius, r ** (-a))  
## Will have a cut off at a minimum size - C-H bond length
#
#min_size = 0.541e-9  # C-H bond length
#max_size = 10000000
##
alpha = 3.0/2  # Possible values: -3/2, -7/2, -19/12 etc. - Mahapatra only say 3/2 and 19/12 produce acceptable ERE spectra.
#norm_const = ((max_size * 1.0e-9) ** (1-alpha) - (1e-9) ** (1-alpha))/(1-alpha)
#sizes = [i * 1.0e-10 for i in list(range(1, max_size))]  # Sizes in 10^{-1} nanometres. 
#size_plot = [math.log(i, 10) for i in sizes]
#numbers = [i ** (-alpha) / norm_const for i in sizes]
#number_plot = [0 if i < 1.0 else math.log(i, 10) for i in numbers]
#
#alpha1 = 5.0/2  # Mahapatra et. al. 2014 dismiss this value as giving realistic spectra.
#norm_const1 = ((max_size * 1.0e-9) ** (1-alpha1) - (1e-9) ** (1-alpha1))/(1-alpha1)
#numbers1 = [i ** (-alpha1) / norm_const1 for i in sizes]
#number_plot1 = [0 if i < 1.0 else  math.log(i, 10) for i in numbers1]
#
#alpha2 = 3.0/2
#norm_const2 = ((max_size * 1.0e-9) ** (1-alpha2) - (1e-9) ** (1-alpha2))/(1-alpha2)
#numbers2 = [i ** (-alpha2) / norm_const2 for i in sizes]
#number_plot2 = [0 if i < 1.0 else  math.log(i, 10) for i in numbers2]
#
#alpha3 = 7.0/2
#norm_const3 = ((max_size * 1.0e-9) ** (1-alpha3) - (1e-9) ** (1-alpha3))/(1-alpha3)
#numbers3 = [i ** (-alpha3) / norm_const3 for i in sizes]
#number_plot3 = [0 if i < 1.0 else math.log(i, 10) for i in numbers3]

#Plotting the size distribution in terms of mass
con = 1.66054e-27
masses = list(range(16, 129)) # in amu
masskg = [math.log(i * con,10) for i in masses] # in kg
#part = [math.log(i ** (- 19.0/12), 10) for i in masses]
#part2 = [math.log(i ** (-1.5),10) for i in masses]
#part4 = [math.log(i ** (-3.5),10) for i in masses]
part = [i ** (- 19.0/12) for i in masses]
part2 = [i ** (-1.5) for i in masses]
part4 = [i ** (-3.5) for i in masses]
norm_mass1 = abs(trapz(part, masses, dx=1))
norm_mass2 = abs(trapz(part2, masses, dx=1))
norm_mass4 = abs(trapz(part4, masses, dx=1))
numberplot = [math.log(i / norm_mass1,10) for i in part]
numberplot2 = [math.log(i / norm_mass2,10) for i in part2]
numberplot4 = [math.log(i / norm_mass4,10) for i in part4]
#
## Using mass of molecules instead of size - # of nucleons * atomic mass unit. 
meth_mass = 16 # all in amu
eth_mass = 30
prop_mass = 44
pent_mass = 72
benz_mass = 78
naph_mass = 128

plt.figure(2)
#pylab.plot(size_plot, number_plot, 'k-')
#pylab.plot(size_plot, number_plot1, 'b-')
#pylab.plot(size_plot, number_plot2, 'g-')
#pylab.plot(size_plot, number_plot3, 'r-')
pylab.plot(masskg, numberplot, 'k-')
#pylab.xlabel("Size (log_10 (m (amu)))")
#pylab.ylabel("Number (log_10 (N))")
#pylab.title("Decay constant: 19/12")
#plt.figure(7)
pylab.plot(masskg,numberplot2,'b-')
#pylab.xlabel("Size (log_10 (m (amu)))")
#pylab.ylabel("Number (log_10 (N))")
#pylab.title("Decay constant: 3/2")
#plt.figure(8)
#pylab.plot(masskg,numbersM3,'k-')
#pylab.xlabel("Size (log_10 (m (amu)))")
#pylab.ylabel("Number (log_10 (N))")
#pylab.title("Decay constant: 5/2")
#plt.figure(9)
pylab.plot(masskg,numberplot4,'g-')
pylab.xlabel("Size (log_10 (m (kg)))")
pylab.ylabel("Number (log_10 (N))")
pylab.title("Size distributions.")

# Effective Z_eff plot
methane = open("methane.txt", "r")
ethane = open("ethane.txt", "r")
propane = csv.reader(open("propane.csv", "r"), delimiter='\t')
pentane = csv.reader(open("pentane.csv", "r"), delimiter='\t')
benzene = csv.reader(open("benzene.csv", "r"), delimiter='\t')
naphthalene = csv.reader(open("naphthalene.csv", "r"), delimiter='\t')
areas = open("areas.txt", "r")
areaLines = [float(line[0:-1]) for line in areas]

methane_levels = [line[0:4] for line in methane]
methane_lvls = [int(i) / 8065.0 for i in methane_levels[4:len(methane_levels)]]
ethane_levels = [line[0:4] for line in ethane]
ethane_lvls = [int(i) / 8065.0 for i in ethane_levels[4:len(ethane_levels)]]
propane_lines = [line for line in propane]
propane_levels = [int(i[2]) / 8065.0 for i in propane_lines]
pentane_lines = [line for line in pentane]
pentane_levels = [int(i[2]) / 8065.0 for i in pentane_lines]
benzene_lines = [line for line in benzene]
benzene_levels = [int(i[2]) / 8065.0 for i in benzene_lines]
naphthalene_lines = [line for line in naphthalene]
naphthalene_levels = [int(i[2]) / 8065.0 for i in naphthalene_lines]
meth_areas = areaLines[0:3]
eth_areas = areaLines[4:15]
prop_areas = areaLines[16:43]
pent_areas = areaLines[44:88]
benz_areas = areaLines[89:108]
naph_areas = areaLines[109:156]

k_B = 8.6173303 * 1e-5  # Boltzmann constant in eV K-1
roomT = 300  # room temperature
HIM_T = 1e6  # hot ionised medium
WIM_T = 8000 # warm ionised medium
WNM_T_min = 6000 # warm neutral medium
WNM_T_max = 10000  # 6000 - 10000 K
CNM_T_min = 50   # cold neutral medium
CNM_T_max = 100    # 50 - 100 K
MM_T_min = 10    # molecular clouds
MM_T_max = 20      # 10 - 20 K
width = 1.5 * k_B * 300 # width of Gaussian
methane_EB = 0  # binding energies
ethane_EB = 0
propane_EB = 10.0e-3
pentane_EB = 60.0 * 1e-3
benzene_EB = 150.0 * 1e-3
naphthalene_EB = 300.0 * 1e-3  # approximately - for some modes the peak is below zero because EB > E_V. 

# Z_eff values
methane_Zeff = 142
ethane_Zeff = 1780
propane_Zeff = 3500
pentane_Zeff = 40200
benzene_Zeff = 15000
naphthalene_Zeff = 494000

### Plotting the energy distribution.
const = 0.3 * 10 ** (-3)
energies = [i / 10.0 for i in range(69)]
Edist_values = [const * i for i in energies]
#Edist = [energy_dist(1000, 10 ** (-5), i, 3.0) for i in energies] 

# To run it like a Monte Carlo - adjusting for probabilities of different energies (but not really as we are not taking random choices of positrons)
# Does provide a good base for doing an actual Monte Carlo in future however.
energies2 = [i * 1.0e-3 for i in range(0, 6800)]# range of positron incident energies, greater than 0, up to 1 eV. 
Edist2 = [const * i for i in energies2]
energy_prob = [int(1000000 * i) for i in Edist2]
j = 0
#weighting = []
#while j < len(Edist2):
#    weighting.extend(energies2[j:j+1] * energy_prob[j]) 
#    j += 1
#print("Done")


def gauss(Zeff, width, E_V, EB):  # E_V is mode energy
    const1 = Zeff / math.sqrt(2 * math.pi * (width ** 2))
    const2 = [((i - E_V + EB) ** 2) / (2 * (width ** 2)) for i in energies2]     # left downshift
    #const2 = [((i - EB) ** 2) / (2 * (width ** 2)) for i in energies2]    # right downshift# with mode energies and left downshift.
    #const2 = [((i - E_V + EB) ** 2) / (2 * (width ** 2)) for i in weighting]  # with energies adjusted for probability
    func = [const1 * math.exp(-i) for i in const2]
    #func = []
    #for i in const2:
    #    if  const2.index(i) <= 6801:
    #        func.append(const1 * math.exp(- i))
    #    else:
    #        func.append(0)
    return func

def listsum(list1, list2):
    '''Performs element-wise addition of lists.'''
    i = 0
    result = []
    while i < len(list1):
        result.append(list1[i] + list2[i])
        i += 1
    return result

def listmult(list1, list2):
    i = 0
    result = []
    while i < len(list1):
        result.append(list1[i] * list2[i])
        i += 1
    return result

alphatest = 1.5 # 1.5, 3.5, 19.0/12
u = 1.66054e-27
# Using mass of nuclei of molecules - ignoring electrons.
#methane_numb = math.log((16 * u) ** (-alphatest))
#ethane_numb = math.log((29 * u) ** (-alphatest))
#propane_numb = math.log((44 * u) ** (-alphatest))
#pentane_numb = math.log((72 * u) ** (-alphatest))
#benzene_numb = math.log((78 * u) ** (-alphatest))
#naph_numb = math.log((128 * u) ** (-alphatest))
methane_numb = (16 * u) ** (-alphatest)
ethane_numb = (29 * u) ** (-alphatest)
propane_numb = (44 * u) ** (-alphatest)
pentane_numb = (72 * u) ** (-alphatest)
benzene_numb = (78 * u) ** (-alphatest)
naph_numb = (128 * u) ** (-alphatest)


# Best to plot each separately.
#methane_sum = [i / meth_areas[0] for i in gauss(methane_Zeff, width, methane_lvls[0], methane_EB)]
#ethane_sum = [i / eth_areas[0] for i in gauss(ethane_Zeff, width, ethane_lvls[0], ethane_EB)]
#propane_sum = [i / prop_areas[0] for i in gauss(propane_Zeff, width, propane_levels[0], propane_EB)]
#pentane_sum = [i / pent_areas[0] for i in gauss(pentane_Zeff, width, pentane_levels[0], pentane_EB)]
#benzene_sum = [i / benz_areas[0] for i in gauss(benzene_Zeff, width, benzene_levels[0], benzene_EB)]
#naph_sum = [i / naph_areas[0] for i in gauss(naphthalene_Zeff, width, naphthalene_levels[0], naphthalene_EB)]
methane_sum = gauss(methane_Zeff, width, methane_lvls[0], methane_EB)
ethane_sum = gauss(ethane_Zeff, width, ethane_lvls[0], ethane_EB)
propane_sum = gauss(propane_Zeff, width, propane_levels[0], propane_EB)
pentane_sum = gauss(pentane_Zeff, width, pentane_levels[0], pentane_EB)
benzene_sum = gauss(benzene_Zeff, width, benzene_levels[0], benzene_EB)
naph_sum = gauss(naphthalene_Zeff, width, naphthalene_levels[0], naphthalene_EB)

# Energy distribution bits
const = 3e-3
energies = [i / 1000.0 for i in range(1,6801)]
Edist_values = [i * const for i in energies]
Edist = [energy_dist(1000, 10 ** (-5), i, 3.0) for i in energies] 


# Normalisation: method 1 - divide gaussian by # of molecules of the relevant species.
plt.figure(6)
k = 1
while k < len(methane_lvls):
    #listsum(methane_sum, [i / meth_areas[k-1] for i in gauss(methane_Zeff, width, methane_lvls[k], methane_EB)])
    listsum(methane_sum, gauss(methane_Zeff, width, methane_lvls[k], methane_EB))
    k += 1
methane_net = [i / methane_numb for i in methane_sum]
#methane_net = [i / (methane_numb * 5.0857497371073558e-32) for i in methane_sum]
#pylab.plot(energies2, methane_net, 'k-')
#pylab.plot(energies2, methane_sum, 'g-')
k = 1
while k < len(ethane_lvls):
    #listsum(ethane_sum, [i / eth_areas[k-1] for i in gauss(ethane_Zeff, width, ethane_lvls[k], ethane_EB)])
    #listsum(ethane_sum, gauss(ethane_Zeff, width, ethane_lvls[k], ethane_EB))
    #pylab.plot(energies2, [i / eth_areas[k-1] for i in gauss(ethane_Zeff, width, ethane_lvls[k], ethane_EB)],'k-')
    k += 1
ethane_net = [i / ethane_numb for i in ethane_sum]
#ethane_net = [i / (ethane_numb * 5.0857497371073558e-32) for i in ethane_sum]
#pylab.plot(energies2, ethane_net, 'k-')
#pylab.plot(energies2, ethane_sum, 'b-')
k = 1
while k < len(propane_levels):
    #listsum(propane_sum, [i / prop_areas[k-1] for i in gauss(propane_Zeff, width, propane_levels[k], propane_EB)])
    listsum(propane_sum, gauss(propane_Zeff, width, propane_levels[k], propane_EB))
    k += 1
propane_net = [i / propane_numb for i in propane_sum]
#propane_net = [i / (propane_numb * 5.0857497371073558e-32) for i in propane_sum]
#pylab.plot(energies2, propane_net, 'k-')
#pylab.plot(energies2, propane_sum, 'r-')
k = 1
while k < len(pentane_levels):
    #listsum(pentane_sum, [i / pent_areas[k-1] for i in gauss(pentane_Zeff, width, pentane_levels[k], pentane_EB)])
    listsum(pentane_sum, gauss(pentane_Zeff, width, pentane_levels[k], pentane_EB))
    k += 1
pentane_net = [i / pentane_numb for i in pentane_sum]
#pentane_net = [i / (pentane_numb * 5.0857497371073558e-32) for i in pentane_sum]
#pylab.plot(energies2, pentane_net, 'k-')
#pylab.plot(energies2, pentane_sum, 'g-')
k = 1
while k < len(benzene_levels):
    #listsum(benzene_sum, [i / benz_areas[k-1] for i in gauss(benzene_Zeff, width, benzene_levels[k], benzene_EB)])
    listsum(benzene_sum, gauss(benzene_Zeff, width, benzene_levels[k], benzene_EB))
    k += 1
benzene_net = [i / benzene_numb for i in benzene_sum]
#benzene_net = [i / (benzene_numb * 5.0857497371073558e-32) for i in benzene_sum]
#pylab.plot(energies2, benzene_net, 'k-')
#pylab.plot(energies2, benzene_sum, 'k+')
k = 1
while k < len(naphthalene_levels):
    #listsum(naph_sum, [i / naph_areas[k-1] for i in gauss(naphthalene_Zeff, width, naphthalene_levels[k], naphthalene_EB)])
    listsum(naph_sum, gauss(naphthalene_Zeff, width, naphthalene_levels[k], naphthalene_EB))
    k += 1
#naph_net = [i / naph_numb for i in naph_sum]
naph_net = [i / naph_numb for i in naph_sum]
#pylab.plot(energies2, naph_net, 'b-')
#pylab.plot(energies2, naph_sum, 'b-')
#mega_sum = listsum(methane_sum, listsum(ethane_sum, listsum(propane_sum, listsum(pentane_sum, listsum(benzene_sum, naph_sum)))))
mega = listsum(methane_net, listsum(ethane_net, listsum(propane_net, listsum(pentane_net, listsum(benzene_net, naph_net)))))
area = trapz(mega,dx=0.001)
mega_sum = [i / area for i in mega]
convolve = listmult(energies, mega_sum)
#mega_sum = [(10 ** 5) * i for i in listsum(methane_net, listsum(ethane_net, listsum(propane_net, listsum(pentane_net, listsum(benzene_net, naph_net)))))]
pylab.plot(energies2, mega_sum, 'r-')
#pylab.plot(weighting, mega_sum, 'r+')
#pylab.plot(energies, Edist_values, 'k-')
#pylab.plot(energies, convolve, 'r-')
#pylab.title("Product of the energy distribution and annihilation spectra" + "\n")
#pylab.title("Annihilation spectra at T = 300 K, normalised" + "\n" + "to area of one and corrected for number abundance.")
#pylab.title("Energy distribution and annihilation spectra (T = 10000 K)" + '\n') 
#pylab.title("Positron resonance with a specific vibrational mode" + "\n" + "of naphthalene in the molecular cloud.")
#pylab.title("Annihilation spectra in hot ionised media (T = 10^6 K)" + "\n")
pylab.xlabel("Positron energy (eV)")
#pylab.ylabel("Z_eff and P(E)")
pylab.ylabel("Z_eff")
pylab.title("Annihilation spectrum for the vibrational modes of ethane" + "\n")
pylab.xlim([0,1.0])
#pylab.ylim([1e-4,1.5e-4])

plt.figure(1)
pylab.plot(energies, Edist_values, 'k+')
#pylab.plot(energies2, mega_sum, 'r-')
pylab.xlabel("Positron kinetic energy (eV)")
pylab.ylabel("P(E) and Z_eff")
pylab.title("Energy distribution and annihilation spectra (T = 6000 K)" + "\n")
methane.close()
ethane.close()

pylab.show()

# To do and method ###

# Testing with different Z_eff, E_b and Gamma_v values (1e-2 to 1e-6) (Gribakin 2010 values).

# Plot Z_eff against total incident energy of positrons.
# Plot for different molecules and obtain a continuous spectrum - effective Zeff plot.

# 1. Have an energy distribution for the positron energies, with probability of that energy on vertical axis. (that from Machacek poster)
# 2. Have a size distribution of the molecule(s), with size (radius or mass) on the horizontal and number on the vertical axis. 
# 3. Use the number of molecules for a given size to get a curve of Z_eff w.r.t. e_b, e and N (pick an incident positron energy here)
# 3a. Use the energy that excites the C-H stretch mode in the molecule (E = omega_v - e_b)
# 4. Convolve the energy distribution from step 1 with another function to obtain a mono-energetic energy distribution.
# 5. Test the plot obtained with different Z_eff values, as well as different e_b values and gamma_v values. 
# 6. Repeat 3 for different molecules (values of N)
# 7. Create a plot of Z_eff against energy for many different molecules at the C-H stretch peak.