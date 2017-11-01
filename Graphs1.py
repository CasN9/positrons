# -*- coding: utf-8 -*-
# Caspian Nicholls, u1027945
# Australian National University

import math
import pylab
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import csv


def energy_dist(Z_eff, gamma, E, E_0):
    '''Creates a Gaussian energy distribution. Z_eff is the amplitude, \
    E_0 is the incident energy, E is the positron energy and gamma is the \
    width of the Gaussian.'''
    coeff = Z_eff / (gamma * math.sqrt(2 * math.pi))
    power1 = - (E - E_0)**2
    power2 = 2* (gamma)**2
    f_E = coeff * math.exp(power1/ power2)
    return f_E

# Plotting the size distribution in terms of mass
conversion = 1.66054e-27
masses = list(range(16, 129)) # in amu
masskg = [math.log(i * conversion,10) for i in masses] # in kg
part = [i ** (-19.0/12) for i in masses]
part2 = [i ** (-1.5) for i in masses]
part3 = [i ** (-3.5) for i in masses]
normalised_mass1 = abs(trapz(part, masses, dx=1))
normalised_mass2 = abs(trapz(part2, masses, dx=1))
normalised_mass4 = abs(trapz(part3, masses, dx=1))
numberplot = [math.log(i / normalised_mass1,10) for i in part]
numberplot2 = [math.log(i / normalised_mass2,10) for i in part2]
numberplot3 = [math.log(i / normalised_mass4,10) for i in part3]

# Using mass of molecules instead of size - # of nucleons * atomic mass unit. 
methane_mass = 16 # all in amu
ethane_mass = 30
propane_mass = 44
pentane_mass = 72
benzene_mass = 78
naphthalene_mass = 128

plt.figure(1)
pylab.plot(masskg, numberplot, 'k-')
pylab.xlabel("Size (log_10 (m (kg)))")
pylab.ylabel("Number (log_10 (N))")
pylab.title("Decay constant: 19/12")
plt.figure(2)
pylab.plot(masskg,numberplot2,'b-')
pylab.xlabel("Size (log_10 (m (kg)))")
pylab.ylabel("Number (log_10 (N))")
pylab.title("Decay constant: 3/2")
plt.figure(3)
pylab.plot(masskg,numberplot3,'g-')
pylab.xlabel("Size (log_10 (m (kg)))")
pylab.ylabel("Number (log_10 (N))")
pylab.title("Decay constant: 7/2")

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
methane_areas = areaLines[0:3]
ethane_areas = areaLines[4:15]
propane_areas = areaLines[16:43]
pentane_areas = areaLines[44:88]
benzene_areas = areaLines[89:108]
naphthalene_areas = areaLines[109:156]

k_B = 8.6173303 * 1e-5  # Boltzmann constant in eV K-1
room_T = 300  # room temperature
HIM_T = 1e6  # hot ionised medium temperature
WIM_T = 8000 # warm ionised medium temperature
WNM_T_min = 6000 # warm neutral medium minimum temperature
WNM_T_max = 10000  # warm neutral medium maximum temperature
CNM_T_min = 50   # cold neutral medium minimum temperature
CNM_T_max = 100    # cold neutral medium maximum temperature
MM_T_min = 10    # molecular cloud minimum temperature
MM_T_max = 20      # molecular cloud maximum temperature
temp = room_T
gaussian_width = 1.5 * k_B * temp
methane_EB = 0  # binding energies
ethane_EB = 0
propane_EB = 10.0e-3
pentane_EB = 60.0 * 1e-3
benzene_EB = 150.0 * 1e-3
naphthalene_EB = 300.0 * 1e-3  # approximately

# Z_eff values
methane_Zeff = 142
ethane_Zeff = 1780
propane_Zeff = 3500
pentane_Zeff = 40200
benzene_Zeff = 15000
naphthalene_Zeff = 494000

### Plotting the energy distribution.
energies2 = [i * 1.0e-3 for i in range(0, 6800)]  # range of positron incident energies, greater than 0, up to 1 eV. 
const = 3e-3
energies = [i / 1000.0 for i in range(1,6801)]
Edist_values = [i * const for i in energies]
Edist = [energy_dist(1000, 10 ** (-5), i, 3.0) for i in energies] 

def gauss(Zeff, width, E_V, EB):  # E_V is mode energy
    const1 = Zeff / math.sqrt(2 * math.pi * (width ** 2))
    const2 = [((i - E_V + EB) ** 2) / (2 * (width ** 2)) for i in energies2]     # left downshift
    func = [const1 * math.exp(-i) for i in const2]
    return func

def listsum(list1, list2):
    '''Performs element-wise addition of two lists.'''
    i = 0
    result = []
    while i < len(list1):
        result.append(list1[i] + list2[i])
        i += 1
    return result

def listmult(list1, list2):
    '''Performs element-wise multiplication of two lists.'''
    i = 0
    result = []
    while i < len(list1):
        result.append(list1[i] * list2[i])
        i += 1
    return result

decay_const = 1.5    # Possible values: 1.5, 3.5, 19.0/12
u = 1.66054e-27      # atomic mass unit, in kg.
methane_numb = (16 * u) ** (-decay_const)
ethane_numb = (29 * u) ** (-decay_const)
propane_numb = (44 * u) ** (-decay_const)
pentane_numb = (72 * u) ** (-decay_const)
benzene_numb = (78 * u) ** (-decay_const)
naphthalene_numb = (128 * u) ** (-decay_const)

methane_sum = gauss(methane_Zeff, gaussian_width, methane_lvls[0], methane_EB)
ethane_sum = gauss(ethane_Zeff, gaussian_width, ethane_lvls[0], ethane_EB)
propane_sum = gauss(propane_Zeff, gaussian_width, propane_levels[0], propane_EB)
pentane_sum = gauss(pentane_Zeff, gaussian_width, pentane_levels[0], pentane_EB)
benzene_sum = gauss(benzene_Zeff, gaussian_width, benzene_levels[0], benzene_EB)
naphthalene_sum = gauss(naphthalene_Zeff, gaussian_width, naphthalene_levels[0], naphthalene_EB)


plt.figure(4)
k = 1
while k < len(methane_lvls):
    listsum(methane_sum, gauss(methane_Zeff, gaussian_width, methane_lvls[k], methane_EB))
    k += 1
methane_net = [i / methane_numb for i in methane_sum]
k = 1
while k < len(ethane_lvls):
    listsum(ethane_sum, gauss(ethane_Zeff, gaussian_width, ethane_lvls[k], ethane_EB))
    pylab.plot(energies2, gauss(ethane_Zeff, gaussian_width, ethane_lvls[k], ethane_EB),'k-')
    pylab.xlabel("Positron energy (eV)")
    pylab.ylabel("Z_eff")
    pylab.title("Annihilation spectrum for the vibrational modes of ethane" + "\n" + "(T = " + str(temp) + " K)")
    pylab.xlim([0,1.0])
    k += 1
ethane_net = [i / ethane_numb for i in ethane_sum]
k = 1
while k < len(propane_levels):
    listsum(propane_sum, gauss(propane_Zeff, gaussian_width, propane_levels[k], propane_EB))
    k += 1
propane_net = [i / propane_numb for i in propane_sum]
k = 1
while k < len(pentane_levels):
    listsum(pentane_sum, gauss(pentane_Zeff, gaussian_width, pentane_levels[k], pentane_EB))
    k += 1
pentane_net = [i / pentane_numb for i in pentane_sum]
k = 1
while k < len(benzene_levels):
    listsum(benzene_sum, gauss(benzene_Zeff, gaussian_width, benzene_levels[k], benzene_EB))
    k += 1
benzene_net = [i / benzene_numb for i in benzene_sum]
k = 1
while k < len(naphthalene_levels):
    listsum(naphthalene_sum, gauss(naphthalene_Zeff, gaussian_width, naphthalene_levels[k], naphthalene_EB))
    k += 1
naphthalene_net = [i / naphthalene_numb for i in naphthalene_sum]
mega = listsum(methane_net, listsum(ethane_net, listsum(propane_net, listsum(pentane_net, listsum(benzene_net, naphthalene_net)))))
area = trapz(mega,dx=0.001)
mega_sum = [i / area for i in mega]
convolve = listmult(energies, mega_sum)

plt.figure(5)
pylab.plot(energies2, mega_sum, 'r-')
pylab.plot(energies, Edist_values, 'k-')
pylab.xlabel("Positron energy (eV)")
pylab.ylabel("Z_eff")
pylab.title("Energy distribution and annihilation spectra (T = " + str(temp) + " K)" + '\n') 
pylab.xlim([0,1.0])

plt.figure(6)
pylab.plot(energies, convolve, 'r-')
pylab.title("Product of the energy distribution and annihilation spectra" \
            + "\n" + "(T = " + str(temp) + " K)")
pylab.ylabel("Z_eff and probability, P(E)") # probability as a function of energy.
pylab.xlabel("Positron energy (eV)")
pylab.xlim([0,1.0])

plt.figure(7)
pylab.plot(energies, Edist_values, 'k+')
pylab.xlabel("Positron kinetic energy (eV)")
pylab.ylabel("Probability, P(E)")
pylab.title("Energy distribution and annihilation spectra (T = " + str(temp) + " K)" + "\n")
pylab.show()

methane.close()
ethane.close()
