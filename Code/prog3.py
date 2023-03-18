# Import libraries 
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import os
import mesa_reader as mr
from scipy.optimize import curve_fit

# PATH DIRECTORY
PATH1_wd = '/Users/tranghuynh/Dropbox/LSU/Fall_2022/ASTR4750/project3/1M_star/LOGS_1_wd/'
PATH15_wd = '/Users/tranghuynh/Dropbox/LSU/Fall_2022/ASTR4750/project3/1M_star/LOGS_15_wd/'
PATH1_H = '/Users/tranghuynh/Dropbox/LSU/Fall_2022/ASTR4750/project3/1M_star/LOGS1_cntr_H_frac/'
PATH15_H = '/Users/tranghuynh/Dropbox/LSU/Fall_2022/ASTR4750/project3/1M_star/LOGS_15_cntr_H_frac/'
f1 = 'history.data'
f15 = 'history.data'
profile1 = 'profile3.data'
profile15 = 'profile1.data'
profile1_wd = 'profile100.data'

################################### BODY ##########################################
# Load data from LOGS using mesa_reader
os.chdir(PATH1_wd)
if os.path.exists(PATH1_wd):
    h = mr.MesaData(f1)
    p1_wd = mr.MesaData(profile1_wd)
    
os.chdir(PATH1_H)
if os.path.exists(PATH1_H):
    p1 = mr.MesaData(profile1)

os.chdir(PATH15_wd)
if os.path.exists(PATH15_wd):
    h15 = mr.MesaData(f15)

os.chdir(PATH15_H)
if os.path.exists(PATH15_H):
    p15 = mr.MesaData(profile15)

    

#Import data for 1 solar mass
ages = h.data('star_age')
log_dt = h.data('log_dt')
log_Teff = h.data('log_Teff')
log_L = h.data('log_L')
rho_c = h.data('center_Rho')
T_c = h.data('center_T')
# print(shape(T_c), '1M')
# print('T_c =', T_c)
# print('\n')

# Model for 15 solar mass

# Import data for 15 solar mass
ages15 = h15.data('star_age')
log_dt15 = h15.data('log_dt')
log_Teff15 = h15.data('log_Teff')
log_L15 = h15.data('log_L')
rho_c15 = h15.data('center_Rho')
T_c15 = h15.data('center_T')
# print(shape(T_c15), '15M')
# print('T_c15 = ', T_c15)


# HR Diagram plot for 1M
plt.figure(1)
plt.plot(h.log_Teff, h.log_L, 'b-', label = ' 1 solar mass')
plt.plot(h15.log_Teff, h15.log_L, 'r-', label = ' 15 solar mass')
plt.title('HR Diagram for 1M and 15M models')
plt.xlabel('log Effective Temperature')
plt.ylabel('log Luminosity')
legend = plt.legend(loc='lower right',shadow='True')
plt.minorticks_on()
# plt.show()


# Evolution of central density versus central temperature
plt.figure(2)
plt.plot(T_c, rho_c,'b-', label = ' 1 solar mass' )
plt.title('Evolution of Central Density versus Central Temperature ')
plt.xlabel(' Central Temperature')
plt.ylabel(' Cenral Density')
legend = plt.legend(loc='lower right',shadow='True')
plt.minorticks_on()
# plt.show()

plt.figure(3)
plt.plot(T_c15, rho_c15,'r-', label = ' 15 solar mass' )
plt.title('Evolution of Central Density versus Central Temperature ')
plt.xlabel(' Central Temperature')
plt.ylabel(' Cenral Density')
legend = plt.legend(loc='lower right',shadow='True')
plt.minorticks_on()
# plt.show()

# PART B
# Load profile data 

# Retrieve profile data for 1M for xa_central_lower_limit(1) = 0.6
log_P = p1.logP
# P = 10 ** log_P
log_rho = p1.logRho
# rho = 10 ** log_rho

#Retrieve profile data for 15M for xa_central_lower_limit(1) = 0.6
log_P15 = p15.logP
log_rho15 = p15.logRho
# P15 = 10 ** log_P15
# rho15 = 10 ** log_rho15

# Retrieve profile data for 1M of white draft stage
log_P_wd = p1_wd.logP
log_rho_wd = p1_wd.logRho


## MAKE PLOT FOR TESTING
# plt.figure(4)
# plt.plot(rho15, P15, label = ' 15 solar mass ')
# plt.xlabel(' Density ')
# plt.ylabel(' Pressure ')
# legend = plt.legend(loc='lower right',shadow='True')
# plt.show()

# #Define EOS function
def eos(log_rho, gamma):
    return log_rho*gamma
# y_mod = eos(trunc_rho15, 0.05) #*1e10
# #Perform the fit	(LINEAR)
fit,cov = curve_fit(eos,log_rho,log_P)
fit15,cov15 = curve_fit(eos, log_rho15, log_P15)
fit_wd, cov_wd = curve_fit(eos, log_rho_wd,log_P_wd)
n1 = 1./(fit - 1.)
n15 = 1./(fit15 - 1.)
n1_wd = 1./(fit_wd - 1.)
print("The polytropic index for 1 solar mass is ", n1 , ", where gamma =", fit, "+/-", cov)
print("The polytropic index for 15 solar mass is ", n15, ", where gamma =", fit15, "+/-", cov15)
print("The polytropic index for 1 solar mass in wd stage is ", n1_wd , ", where gamma =", fit_wd, "+/-", cov_wd)


#Fit Plot for 1M with central limit = 0.6
plt.figure(4)
plt.plot(log_rho, log_P, 'r-', label = 'data')
plt.plot(log_rho, eos(log_rho, abs(fit)), 'b-', label = 'fit')
plt.title('Equation of State for 1M when x_central_mass_fraction > 0.6')
plt.xlabel(' Log Density ')
plt.ylabel(' Log Pressure ')
legend = plt.legend(loc='lower right',shadow='True')
# plt.show()

#Fit Plot for 15M with central limit = 0.6
plt.figure(5)
plt.plot(log_rho15, log_P15, 'r-', label = 'data')
plt.plot(log_rho15, eos(log_rho15, abs(fit15)), 'b-', label = 'fit')
plt.title('Equation of State for 15M when x_central_mass_fraction > 0.6')
plt.xlabel(' Log Density ')
plt.ylabel(' Log Pressure ')
legend = plt.legend(loc='lower right',shadow='True')

#Fit Plot for 1M at white dwarf stage
plt.figure(6)
plt.plot(log_rho_wd, log_P_wd, 'r-', label = 'data')
plt.plot(log_rho_wd, eos(log_rho_wd, abs(fit_wd)), 'b-', label = 'fit')
plt.title('Equation of State for 1M at white dwarf stage')
plt.xlabel(' Log Density ')
plt.ylabel(' Log Pressure ')
legend = plt.legend(loc='lower right',shadow='True')
plt.show()








