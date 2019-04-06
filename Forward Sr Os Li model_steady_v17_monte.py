#Forward Sr Os Li model
"""
Created on Thu Jan 31 15:25:22 2019

@author: joach
"""
import xlsxwriter
import scipy as sci 
from scipy.integrate import solve_ivp
import sympy
import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
import scipy.interpolate
import matplotlib.pyplot as plt
import math
import time
import pickle
import dill
from scipy.signal import savgol_filter
import seaborn as sns


#assign symbolic variables
dLi_ocean, der_dLi_ocean, dOs_ocean, der_dOs_ocean, dSr_ocean, der_dSr_ocean, K_riv, K_ht, K_lt, K_dia, K_light,\
K_dust, K_cosm = sympy.symbols("dLi_ocean der_dLi_ocean dOs_ocean der_dOs_ocean dSr_ocean der_dSr_ocean K_riv K_ht K_lt K_dia K_light K_dust K_cosm")


#time at model start and end
AGE_OLD = 201.6
AGE_YOUNG = 199

#osmium molar mass
Os_mol = 190.2 #g/mol

#strontium molar mass
Sr_mol = 87.62 #g/mol

#lithium molar mass
Li_mol = 6.941 #g/mol


#time step in yrs
t_interval = 1000

#number of resampling
num_monte = 20

num_points = int(np.floor((AGE_OLD - AGE_YOUNG) * 1E6 / t_interval) + 1)

#time length of steady state, or time before perturbation
pre_interval = 1E5

#t_pre = [0, pre_interval - t_interval]
#t_pre = [0, 299000]

#time length of model to run
t_end = 2.6E6


def perturb(mag, period):
    period = int(period)
    mag = mag - 1
    x_temp = [0, pre_interval - 1000]
    y_temp = [1,1]
    for x in range(0, period + 1, t_interval):
        y = (mag * np.sin(x/period * np.pi)) + 1
        x = x + pre_interval
        x_temp.append(x)
        y_temp.append(y)
    x_temp.append(t_end)
    y_temp.append(1)
    K_result = y_temp
    t_result = x_temp
    x_temp = []
    y_temp = []
    return t_result, K_result

#Light flux of Li will have to have different distribution since it's K should be 0 before and after perturbation
def light_perturb(mag, period):
    period = int(period)
    mag = mag
    x_temp = [0, pre_interval - 1000]
    y_temp = [0,0]
    for x in range(0, period + 1, t_interval):
        y = (mag * np.sin(x/period * np.pi))
        x = x + pre_interval
        x_temp.append(x)
        y_temp.append(y)
    x_temp.append(t_end)
    y_temp.append(0)
    K_result = y_temp
    t_result = x_temp
    x_temp = []
    y_temp = []
    return t_result, K_result


def Li_lt_perturb(mag, period):
    period = int(period)
    mag = mag - 13
    x_temp = [0, pre_interval - 1000]
    y_temp = [13,13]
    for x in range(0, period + 1, t_interval):
        y = (mag * np.sin(x/period * np.pi)) + 13
        x = x + pre_interval
        x_temp.append(x)
        y_temp.append(y)
    x_temp.append(t_end)
    y_temp.append(13)
    K_result = y_temp
    t_result = x_temp
    x_temp = []
    y_temp = []
    return t_result, K_result

    
def temp_perturb(mag, period):
    period = int(period)
    mag = mag - 4
    x_temp = [0, pre_interval - 1000]
    y_temp = [4,4]
    for x in range(0, period + 1, t_interval):
        y = (mag * np.sin(x/period * np.pi)) + 4
        x = x + pre_interval
        x_temp.append(x)
        y_temp.append(y)
    x_temp.append(t_end)
    y_temp.append(4)
    K_result = y_temp
    t_result = x_temp
    x_temp = []
    y_temp = []
    return t_result, K_result



t_change1 = np.arange(0, t_end, 1000)



#Os mass balance parameters, all fluxes in mol/yr
# Fluxes and per mil from Misra and Froelich, 2012; Pogge et al, 2013; Coogan et al, 2017; Li and Elderfield, 2013

Os_ht = 10.5 #high Temp hydrothermal Os flux
R_Os_ht = 0.26 #high Temp hydrothermal Os ratio

Os_lt = 105.1 #low Temp hydrothermal Os flux
R_Os_lt = 0.11 #low Temp hydrothermal Os ratio

Os_cosm = 52.6 #cosmogenic Os flux
R_Os_cosm = 0.127 #cosmogenic Os ratio

Os_dust = 36.8 #dust Os flux
R_Os_dust = 1.05 #dust Os ratio

Os_riv = 1800 #Riverine Os flux
R_Os_riv = 0.8 #Riverine Os ratio

#normalizing ratios to get more accurate calculation (see Li and Elderfield 2013)
R_Os_ht_norm = R_Os_ht / (7.4 + R_Os_ht)
R_Os_lt_norm = R_Os_lt / (7.4 + R_Os_lt)
R_Os_cosm_norm = R_Os_cosm / (7.4 + R_Os_cosm)
R_Os_dust_norm = R_Os_dust / (7.4 + R_Os_dust)
R_Os_riv_norm = R_Os_riv / (7.4 + R_Os_riv)

#Concentration of Os in ocean in mol
Os_ocean = 7.20e7

# Sr mass balance parameters, all fluxes in mol/yr
# Fluxes and per mil from Misra and Froelich, 2012; Pogge et al, 2013; Coogan et al, 2017; Li and Elderfield, 2013
Sr_dia = 3.4e9 
R_Sr_dia = 0.7084

Sr_ht = 8.4e9
R_Sr_ht = 0.7025

#Need to change!
Sr_lt = 26e9
k_ini = 10E15* math.exp(-92/((8.314472/1000) * 286.15)) #from Coogan and Dosso, 2015 EPSL using initial bottom water temp of 4 deg C plus 9 deg C (see paper) and converted to Kelvin = 286.15
R_Sr_lt = math.exp(-k_ini)*(dSr_ocean - 0.7025) + 0.7025

Sr_riv = 33.7e9
R_Sr_riv = 0.71144 #value from Vance et al 2009

#normalizing ratios to get more accurate calculation (see Li and Elderfield 2013)
R_Sr_dia_norm = R_Sr_dia / (R_Sr_dia + 9.43)
R_Sr_ht_norm = R_Sr_ht / (R_Sr_ht + 9.43)
R_Sr_riv_norm = R_Sr_riv / (R_Sr_riv + 9.43)

#Concentration of Sr in ocean in mol
Sr_ocean = 1.25e17

# Li mass balance parameters, all fluxes in mol/yr, ratios in permil
# Fluxes and per mil from Misra and Froelich, 2012; Pogge et al, 2013; Coogan et al, 2017; Li and Elderfield, 2013
Li_riv = 10e9 
R_Li_riv = 23

# note: sz the light of hydrothermal is different than he light of rocks
Li_ht = 5.2e9
R_Li_ht = 6.3

Li_lt = 15.2e9
R_Li_lt = 13


Li_light = 10e9
R_Li_light = 1

#normalizing ratios to get more accurate calculation (see Li and Elderfield 2013)
R_Li_riv_norm = R_Li_riv / (R_Li_riv + 1)
R_Li_ht_norm = R_Li_ht / (R_Li_ht + 1)
R_Li_lt_norm = R_Li_lt / (R_Li_lt + 1)
R_Li_light_norm = R_Li_light / (R_Li_light + 1)

#Multiplied conc of Li in seawater (0.000026 M) times V of seawater from Li and Elderfield 2013
Li_ocean = 3.6E16

#Li_ocean = 789856689309622 #Reduced N in this version by multiplying residence time (2 Myr) * riverine Li flux


#mass balance equations for Sr, Os, and Li at steady-state

f_dSr_norm = sympy.Eq(K_ht * Sr_ht * (R_Sr_ht - dSr_ocean) + K_lt * Sr_lt * (R_Sr_lt - dSr_ocean) + K_dia * Sr_dia * (R_Sr_dia - dSr_ocean) + K_riv * Sr_riv * (R_Sr_riv - dSr_ocean) - Sr_ocean * der_dSr_ocean, 0)

f_dOs_norm = sympy.Eq(K_ht * Os_ht * (R_Os_ht - dOs_ocean) + K_lt * Os_lt * (R_Os_lt - dOs_ocean) + K_cosm * Os_cosm * (R_Os_cosm - dOs_ocean) + K_dust * Os_dust * (R_Os_dust - dOs_ocean) + K_riv * Os_riv * (R_Os_riv - dOs_ocean) - Os_ocean * der_dOs_ocean, 0)
    
f_dLi_norm = sympy.Eq(K_ht * Li_ht * (R_Li_ht - dLi_ocean) - K_lt * Li_lt * ((dLi_ocean - R_Li_lt) - dLi_ocean) + K_riv * Li_riv * (R_Li_riv - dLi_ocean) + K_light * Li_light * (R_Li_light - dLi_ocean) - Li_ocean * der_dLi_ocean, 0)


sol = sympy.solve([f_dOs_norm, f_dSr_norm, f_dLi_norm], [K_ht, K_lt, K_riv])

F_K_ht = sympy.simplify(sol[K_ht])
F_K_lt = sympy.simplify(sol[K_lt])
F_K_riv = sympy.simplify(sol[K_riv])

##F_K_light = sympy.simplify(sol[K_light])
##F_K_dust = sympy.simplify(sol[K_dust])
##F_K_dia = sympy.simplify(sol[K_dia])
##F_K_cosm = sympy.simplify(sol[K_cosm])
#
F_K_ht = sympy.lambdify((dOs_ocean, dSr_ocean, dLi_ocean, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm), F_K_ht)
F_K_lt = sympy.lambdify((dOs_ocean, dSr_ocean, dLi_ocean, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm), F_K_lt)
F_K_riv = sympy.lambdify((dOs_ocean, dSr_ocean, dLi_ocean, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm), F_K_riv)

##Creates matices to store parameters
#K_ht_total = np.ones((num_monte, 1)) * -1
#K_lt_total = np.ones((num_monte, 1)) * -1
#K_riv_total = np.ones((num_monte, 1)) * -1

#Load in isotope records to construct ODE
#Later, will use isotope records to constrain outputs:if result at each time step is within a certain percent range away from true value, then keep it. And record success rate

df_Sp = pd.read_excel("TJ_data.xlsx", sheet_name = "Sp")
age_Sp = df_Sp['Age (Ma)']
Sp_record = df_Sp['Sp']

df_TJ_volc = pd.read_excel("TJ_data.xlsx", sheet_name = "Volc")
age_TJ_volc = df_TJ_volc['Age (Ma)']
TJ_volc_record = df_TJ_volc['Volc']

df_dC_carbb = pd.read_excel("TJ_data.xlsx", sheet_name = "dC_carbb")
age_dC_carbb = df_dC_carbb['Age (Ma)']
dC_carbb_record = df_dC_carbb['dC_carbb']

df_dC_orgb = pd.read_excel("TJ_data.xlsx", sheet_name = "dC_orgb")
age_dC_orgb = df_dC_orgb['Age (Ma)']
dC_orgb_record = df_dC_orgb['dC_orgb']

df_dLi_ocean = pd.read_excel("TJ_data.xlsx", sheet_name = "dLi_ocean")
age_dLi_ocean = df_dLi_ocean['Age (Ma)']
dLi_ocean_record = df_dLi_ocean['dLi_ocean']

df_dSr_ocean = pd.read_excel("Italy_TJ_Sr.xlsx", sheet_name = "Sheet3")
age_dSr_ocean = df_dSr_ocean['Age (Ma)']
dSr_ocean_record = df_dSr_ocean['dSr_ocean']

#df_dSr_ocean = pd.read_excel("TJ_data.xlsx", sheet_name = "dSr_ocean")
#age_dSr_ocean = df_dSr_ocean['Age (Ma)']
#dSr_ocean_record = df_dSr_ocean['dSr_ocean']

df_dOs_ocean = pd.read_excel("TJ_data.xlsx", sheet_name = "dOs_ocean")
age_dOs_ocean = df_dOs_ocean['Age (Ma)']
dOs_ocean_record = df_dOs_ocean['dOs_ocean']

df_co2 = pd.read_excel("TJ_data.xlsx", sheet_name = "CO2")
age_co2 = df_co2['Age (Ma)']
co2_record = df_co2['CO2_mean (ppm)']
co2_record_yerr_low = df_co2['yerr_low']
co2_record_yerr_high = df_co2['yerr_high']

df_T = pd.read_excel("TJ_data.xlsx", sheet_name = "T")
age_T = df_T['Age (Ma)']
T_record = df_T['T_mean']
T_record_low = df_T['T_min']
T_record_high = df_T['T_max']
T_record_yerr_low = T_record -T_record_low
T_record_yerr_high = T_record_high - T_record

#interpolating the data record

#for spreading rate of seafloor
age_Sp = (AGE_OLD - age_Sp) * 1E6
Sp_interp = sci.interpolate.interp1d(age_Sp, Sp_record, kind = "cubic")

#for lithium isotopes
age_dLi_ocean = (AGE_OLD - age_dLi_ocean) * 1E6
dLi_ocean_record_smooth = lowess(dLi_ocean_record, age_dLi_ocean, missing = "drop")
dLi_ocean_interp = sci.interpolate.interp1d(dLi_ocean_record_smooth[:,0], dLi_ocean_record_smooth[:,1], kind = "cubic")

#For carbon isotopes in carbonates
age_dC_carbb = (AGE_OLD - age_dC_carbb) * 1E6
dC_carbb_record_smooth = lowess(dC_carbb_record, age_dC_carbb, missing = "drop")
dC_carbb_interp = sci.interpolate.interp1d(dC_carbb_record_smooth[:,0], dC_carbb_record_smooth[:,1], kind = "cubic")

#For carbon isotopes in org C
age_dC_orgb = (AGE_OLD - age_dC_orgb)* 1E6
dC_orgb_record_smooth = lowess(dC_orgb_record, age_dC_orgb, missing = "drop")
dC_orgb_interp = sci.interpolate.interp1d(dC_orgb_record_smooth[:,0], dC_orgb_record_smooth[:,1], kind = "cubic")

#For Sr isotopes in seawater
age_dSr_ocean = (AGE_OLD - age_dSr_ocean) * 1E6
dSr_ocean_record_smooth = lowess(dSr_ocean_record, age_dSr_ocean, missing = "drop")
dSr_ocean_interp = sci.interpolate.interp1d(dSr_ocean_record_smooth[:,0], dSr_ocean_record_smooth[:,1], kind = "cubic")

#For Os isotopes in seawater
age_dOs_ocean = (AGE_OLD - age_dOs_ocean) * 1E6
dOs_ocean_record_smooth = lowess(dOs_ocean_record, age_dOs_ocean, missing = "drop")
dOs_ocean_interp = sci.interpolate.interp1d(dOs_ocean_record_smooth[:,0], dOs_ocean_record_smooth[:,1], kind = "cubic")

#For CO2 record
age_co2 = (AGE_OLD - age_co2) * 1E6
co2_record_smooth = lowess(co2_record, age_co2, missing = "drop")
co2_interp = sci.interpolate.interp1d(co2_record_smooth[:,0], co2_record_smooth[:,1], kind = "cubic")

co2_record_yerr_low_smooth = lowess(co2_record_yerr_low, age_co2, missing = "drop")
co2_yerr_low_interp = sci.interpolate.interp1d(co2_record_yerr_low_smooth[:,0], co2_record_yerr_low_smooth[:,1], kind = "cubic")

co2_record_yerr_high_smooth = lowess(co2_record_yerr_high, age_co2, missing = "drop")
co2_yerr_high_interp = sci.interpolate.interp1d(co2_record_yerr_high_smooth[:,0], co2_record_yerr_high_smooth[:,1], kind = "cubic")

#For T record
age_T = (AGE_OLD - age_T) * 1E6
T_record_smooth = lowess(T_record, age_T, missing = "drop")
# T_interp = sci.interpolate.interp1d(age_T, T_record_smooth, kind = "cubic")
T_interp = sci.interpolate.interp1d(T_record_smooth[:,0], T_record_smooth[:,1], kind = "cubic")

T_record_yerr_low_smooth = lowess(T_record_yerr_low, age_T, missing = "drop")
T_yerr_low_interp = sci.interpolate.interp1d(T_record_yerr_low_smooth[:,0], T_record_yerr_low_smooth[:,1], kind = "cubic")

T_record_yerr_high_smooth = lowess(T_record_yerr_high, age_T, missing = "drop")
T_yerr_high_interp = sci.interpolate.interp1d(T_record_yerr_high_smooth[:,0], T_record_yerr_high_smooth[:,1], kind = "cubic")


#1 sigma errors
e_dOs_ocean = 0.05
e_dLi_ht = 0.1 
e_Li_ocean = 0.2  #Multiplied conc of Li in seawater times V of seawater from Li and Elderfield 2013
e_Li_riv = 0.05
e_dLi_ocean = 0.1
e_flux_light_Li = 0.1 

#1 sigma
#For figure 1
e_dC_carbb = 0.2
e_dC_orgb = 0.5
e_dSr_ocean = 0.00002
e_dOs_ocean = 0.05
e_Sp = 0.04

#Intial K values, or changes in the relative intensity of each flux term (both source and sink)
# K = 1 means no change in flux
K_riv = 1 #riverine intensity
K_ht = 1 #high temp hydrothermal intensity
K_lt = 1 #low temp hydrothermal intensity/marine weathering
K_dia = 1 #diagenetic intensity
K_light = 0 #light Li intensity, only turns on (i.e. to 1) during C injection
K_dust = 1 #dust intensity
K_cosm = 1 #cosmogenic intensity


#Model initialization, i.e. solve fluxes at time zero 
t = 0

Sp_now = Sp_interp(t)

dLi_ocean_now = dLi_ocean_interp(dLi_ocean_record_smooth[0,0])

dC_carbb_now = dC_carbb_interp(t)

dC_orgb_now = dC_orgb_interp(t)

dOs_ocean_now = dOs_ocean_interp(dOs_ocean_record_smooth[0,0])

dSr_ocean_now = dSr_ocean_interp(t)

dSr_ocean_old = dSr_ocean_interp(t - t_interval)

dLi_ocean_old = dLi_ocean_interp(dLi_ocean_record_smooth[0,0])

dOs_ocean_old = dOs_ocean_interp(dOs_ocean_record_smooth[0,0])

 
#der_dSr_ocean_monte = (dSr_ocean_now_monte - dSr_ocean_old_monte)/ t_interval
#der_dOs_ocean_monte = (dOs_ocean_now_monte - dOs_ocean_old_monte)/ t_interval
#der_dLi_ocean_monte = (dLi_ocean_now_monte - dLi_ocean_old_monte)/ t_interval
#


# derivative of dOs in ocean is set to zero
der_dOs_ocean = 0
# derivative of dLi in ocean is set to zero
der_dLi_ocean = 0
der_dSr_ocean = 0

##K values for fluxes at time zero (i.e. before perturbation)
#K_ht_total = F_K_ht(dOs_ocean_now_monte, dSr_ocean_now_monte, dLi_ocean_now_monte, der_dOs_ocean_monte, der_dSr_ocean_monte, der_dLi_ocean_monte, K_light, K_dust, K_dia, K_cosm)
#    
#K_lt_total = F_K_lt(dOs_ocean_now_monte, dSr_ocean_now_monte, dLi_ocean_now_monte, der_dOs_ocean_monte, der_dSr_ocean_monte, der_dLi_ocean_monte, K_light, K_dust, K_dia, K_cosm)
#    
#K_riv_total = F_K_riv(dOs_ocean_now_monte, dSr_ocean_now_monte, dLi_ocean_now_monte, der_dOs_ocean_monte, der_dSr_ocean_monte, der_dLi_ocean_monte, K_light, K_dust, K_dia, K_cosm)


#K values for fluxes at time zero (i.e. before perturbation)
K_ht_total = F_K_ht(dOs_ocean_now, dSr_ocean_now, dLi_ocean_now, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm)
    
K_lt_total = F_K_lt(dOs_ocean_now, dSr_ocean_now, dLi_ocean_now, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm)
    
K_riv_total = F_K_riv(dOs_ocean_now, dSr_ocean_now, dLi_ocean_now, der_dOs_ocean, der_dSr_ocean, der_dLi_ocean, K_light, K_dust, K_dia, K_cosm)


print("K_ht_total is: ", K_ht_total, "K_lt_total is: ", K_lt_total, "K_riv_total", K_riv_total)


#Replace any negative values with nan, and only take average of positive values
#K_ht_total[K_ht_total < 0] = np.nan
K_ht_ave = np.nanmean(K_ht_total)

#K_lt_total[K_lt_total < 0] = np.nan
K_lt_ave = np.nanmean(K_lt_total)

#K_riv_total[K_riv_total < 0] = np.nan
K_riv_ave = np.nanmean(K_riv_total)

#Old fluxes calculated via inverse solution
Old_Os_ht = K_ht_total * Os_ht
Old_Os_lt = K_lt_total * Os_lt
Old_Os_riv = K_riv_total * Os_riv

Old_Sr_ht = K_ht_total * Sr_ht 
Old_Sr_lt = K_lt_total * Sr_lt
Old_Sr_riv = K_riv_total * Sr_riv

Old_Li_ht = K_ht_total * Li_ht 
Old_Li_lt = K_lt_total * Li_lt
Old_Li_riv = K_riv_total * Li_riv

Old_Li_light = K_riv_total * Li_light

Os_ocean = 10000 * Old_Os_riv
Sr_ocean = 3E6 * Old_Sr_riv
Li_ocean = 2E6 *Old_Li_riv

t=0

#Draws intensity fluxes from a uniform distribution

#K_dia_monte = np.random.uniform(0, 1.5, num_monte)
#K_light_monte = np.random.uniform(0., 1., num_monte) 
#K_dust_monte = np.random.uniform(0, 1.5, num_monte) 
#K_cosm_monte = np.random.uniform(0, 1.5, num_monte)
#K_riv_monte = 1
#K_ht_monte = 1
#K_lt_monte = 5
#K_dia_monte = 1
#K_light_monte = 0
#K_dust_monte = 1
#K_cosm_monte = 1


#time length of perturbation
perturb_interval = 4E5

K_riv_mag = 10
K_ht_mag = 2
K_lt_mag = 10
K_light_mag = 10
#
#K_riv_orig = perturb(K_riv_mag, perturb_interval)
#K_ht_orig = perturb(K_ht_mag, perturb_interval)
#K_lt_orig = perturb(K_lt_mag, perturb_interval)
#K_light_orig = light_perturb(K_light_mag, perturb_interval)
#
#K_riv_interp = sci.interpolate.interp1d(K_riv_orig[0], K_riv_orig[1])
#K_ht_interp = sci.interpolate.interp1d(K_ht_orig[0], K_ht_orig[1])
#K_lt_interp = sci.interpolate.interp1d(K_lt_orig[0], K_lt_orig[1])
#K_light_interp = sci.interpolate.interp1d(K_light_orig[0], K_light_orig[1])

R_Li_lt_orig = Li_lt_perturb(10, perturb_interval)

#change in bottom water temperature
temp_change = temp_perturb(13, perturb_interval)

#camp = light_perturb(1E7, 6E5)
#camp = light_perturb(5E6, 4E6)
camp = light_perturb(5E6, 1E6) #CAMP activity may have lasted 600 kyr according to Blackburn et al. 2013 Science
f_sw = 115*1000 #115 t/km^2 *yr value from Cohen and Coe P^3 2007 converted to kg/km^2 * yr

R_Li_lt_interp = sci.interpolate.interp1d(R_Li_lt_orig[0], R_Li_lt_orig[1])

camp_interp = sci.interpolate.interp1d(camp[0], camp[1])

temp_interp = sci.interpolate.interp1d(temp_change[0], temp_change[1])


Sr_camp = (f_sw * 1E-3 * 203.7)/Sr_mol #203.7 ppm (mg/kg) Sr is average Sr conce from CAMP basalts from Deckart et al 2005
Os_camp = (f_sw * 1E-3 * 0.00007)/Os_mol #0.00007 is ppm(mg/kg), but Cohen and Coe report 70ppt is the Os conc of CAMP basalts
Li_camp = (f_sw * 1E-3 * 10.7)/Li_mol

R_Sr_camp = 0.7025
R_Os_camp = 0.12
R_Li_camp = 5

fluid_temp = temp_change[1]    
kelvin_temp = [x + 9 + 273.15 for x in fluid_temp] # add 9 because fluid in off-axis system is ~ 9 deg C warmer than bottom seawater (Coogan and Dosso 2015) and then convert from deg C to kelvin
k = [10E15* math.exp(-92/((8.314472/1000) * x)) for x in kelvin_temp]
frac_bas = [math.exp(x) - 1 for x in k]
frac_bas_interp = sci.interpolate.interp1d(temp_change[0], frac_bas)

K_riv_monte = np.random.uniform(1, K_riv_mag, num_monte) 
K_ht_monte = np.random.uniform(1, K_ht_mag, num_monte) 
K_lt_monte = np.random.uniform(1, K_lt_mag, num_monte) 


age_all = np.linspace((-1E5+1), 1.09E6, 2600)
age_all_plot =  AGE_OLD - age_all/10**6
dOs_plot_real = dOs_ocean_interp(age_all)
dLi_plot_real = dLi_ocean_interp(age_all)
dSr_plot_real = dSr_ocean_interp(age_all)

success = pd.DataFrame()
failure = pd.DataFrame()


for i in range(0, num_monte):
    K_riv_orig = perturb(K_riv_monte[i], perturb_interval)
    K_ht_orig = perturb(K_ht_monte[i], perturb_interval)
    K_lt_orig = perturb(K_lt_monte[i], perturb_interval)
    K_light_orig = light_perturb(K_light_mag, perturb_interval)
    
    K_riv_interp = sci.interpolate.interp1d(K_riv_orig[0], K_riv_orig[1])
    K_ht_interp = sci.interpolate.interp1d(K_ht_orig[0], K_ht_orig[1])
    K_lt_interp = sci.interpolate.interp1d(K_lt_orig[0], K_lt_orig[1])
    K_light_interp = sci.interpolate.interp1d(K_light_orig[0], K_light_orig[1])
    
    #Os set up
    K_Os_flux_times_ratio = lambda t, y: K_ht_interp(t)* Old_Os_ht * (R_Os_ht - y) + K_lt_interp(t) * Old_Os_lt * (R_Os_lt - y) + Os_cosm * (R_Os_cosm - y) + Os_dust * (R_Os_dust - y) + K_riv_interp(t) * Old_Os_riv * (R_Os_riv - y) + camp_interp(t) * Os_camp * (R_Os_camp - y)
        
    K_Os_slope_ratio = lambda t, y: K_Os_flux_times_ratio(t,y)/Os_ocean
        
    K_Os_ode_result = solve_ivp(K_Os_slope_ratio, (0, t_end), dOs_ocean_now.flatten(), method = 'LSODA', t_eval = t_change1)
    
    dOs_total = K_Os_ode_result.y
    
  
    dOs_plot = dOs_total.flatten()
    
    #Li set up
    K_Li_flux_times_ratio = lambda t, y: K_ht_interp(t) * Old_Li_ht * (R_Li_ht - y) - K_lt_interp(t) * Old_Li_lt * ((y - R_Li_lt_interp(t)) - y) + K_riv_interp(t) * Old_Li_riv * (R_Li_riv - y) + K_light_interp(t) * Old_Li_light * (R_Li_light - y) + camp_interp(t) * Li_camp * (R_Li_camp - y)
        
    K_Li_slope_ratio = lambda t, y: K_Li_flux_times_ratio(t,y)/Li_ocean
        
    K_Li_ode_result = solve_ivp(K_Li_slope_ratio, (0, t_end) , dLi_ocean_now.flatten(), method = 'LSODA', t_eval = t_change1)
    
    dLi_total = K_Li_ode_result.y
    
    dLi_plot = dLi_total.flatten()
    
    
    #Sr set up
    K_Sr_flux_times_ratio = lambda t, y: K_ht_interp(t) * Old_Sr_ht * (R_Sr_ht - y) + K_lt_interp(t) * Old_Sr_lt * ((1 * y + frac_bas_interp(t) * (0.7025))/(1 + frac_bas_interp(t)) - y) + Sr_dia * (R_Sr_dia - y) + K_riv_interp(t) * Old_Sr_riv * (R_Sr_riv - y) + camp_interp(t) * Sr_camp * (R_Sr_camp - y)
        
    K_Sr_slope_ratio = lambda t, y: K_Sr_flux_times_ratio(t,y)/Sr_ocean
        
    K_Sr_ode_result = solve_ivp(K_Sr_slope_ratio, (0, t_end), dSr_ocean_now.flatten(), method = 'LSODA', t_eval = t_change1)
    
    dSr_total = K_Sr_ode_result.y
    
    dSr_plot = dSr_total.flatten()
    
    
    if all(dLi_plot <= dLi_plot_real + 0.33*dLi_plot_real) and all(dLi_plot >= dLi_plot_real - 0.33*dLi_plot_real):
        data = pd.DataFrame({'K_lt': [K_lt_monte[i]], 'K_ht': [K_ht_monte[i]], 'K_riv': [K_riv_monte[i]], 'K_light': [K_light_mag]})
        success = success.append(data)
    else:
        data_fail = pd.DataFrame({'K_lt': [K_lt_monte[i]], 'K_ht': [K_ht_monte[i]], 'K_riv': [K_riv_monte[i]], 'K_light': [K_light_mag]})
        failure = failure.append(data_fail)


t_plot = t_change1
##This test is to make the actual Os record start from zero
##test = dOs_ocean_record_smooth[:,0] + abs(dOs_ocean_record_smooth[0,0])
#age_dOs_plot = df_dOs_ocean['Age (Ma)'].sort_values(ascending=False)
#age_dOs_plot_real = ((age_dOs_plot[0] - age_dOs_plot)* 1E6)
#
#
#age_dLi_plot = df_dLi_ocean['Age (Ma)'].sort_values(ascending=False)
#age_dLi_plot_real = ((age_dLi_plot[18] - age_dLi_plot)* 1E6)
#
#age_dSr_plot = df_dSr_ocean['Age (Ma)'].sort_values(ascending=False)
#age_dSr_plot_real = ((age_dSr_plot[253] - age_dSr_plot)* 1E6)

upper_dLi_plot_real = dLi_plot_real + 0.33*dLi_plot_real
lower_dLi_plot_real = dLi_plot_real - 0.33*dLi_plot_real

upper_dSr_plot_real = dSr_plot_real + 0.0005*dSr_plot_real
lower_dSr_plot_real = dSr_plot_real - 0.0005*dSr_plot_real

upper_dOs_plot_real = dOs_plot_real + 0.1*dOs_plot_real
lower_dOs_plot_real = dOs_plot_real - 0.1*dOs_plot_real


sns.set_color_codes()
sns.set_style('white')

f, axarr = plt.subplots(nrows = 3, ncols = 2)
f.suptitle(str(int(perturb_interval)) + ' kyr perturbation w/ ' + 'K_lt: ' + str(K_lt_mag) + ', K_ht: ' + str(K_ht_mag) + ', K_riv: ' + str(K_riv_mag) + ', K_light: ' + str(K_light_mag))
axarr[0,0].plot(age_all_plot, dOs_plot, 'r--', label = 'dOs')
#axarr[0,0].plot(t_plot, dOs_plot, 'r--', label = 'dOs')
#axarr[0].set_ylim(0.4,1)
#axarr[0,0].set_xlim(0, 2.6E6)
axarr[0,0].invert_xaxis()
axarr[0,0].legend(loc = 'best')
axarr[0,0].set_title('Model Data')
axarr[0,0].set_xticks([])

axarr[1,0].plot(age_all_plot, dLi_plot, 'g--', label = 'dLi')
#axarr[1,0].plot(t_plot, dLi_plot, 'g--', label = 'dLi')
#axarr[1,0].set_xlim(0, 2.6E6)
axarr[1,0].invert_xaxis()
#axarr[1].set_ylim(10,20)
axarr[1,0].legend(loc = 'best')
axarr[1,0].set_xticks([])

axarr[2,0].plot(age_all_plot, dSr_plot, 'b--', label = 'dSr')
#axarr[2,0].plot(t_plot, dSr_plot, 'b--', label = 'dSr')
#axarr[2,0].set_xlim(0, 2.6E6)
axarr[2,0].invert_xaxis()
#axarr[2].set_ylim(0.7073,0.7078)
axarr[2,0].legend(loc = 'best')

#Real Data
axarr[0,1].plot(age_all_plot, dOs_plot_real, 'r--', label = 'dOs')
#
axarr[0,1].invert_xaxis()
axarr[0,1].legend(loc = 'best')
axarr[0,1].set_title('Real Data')
axarr[0,1].fill_between(age_all_plot, upper_dOs_plot_real, lower_dOs_plot_real, facecolor='red', alpha=0.3)
#axarr[0,1].set_xlim(0, 5e6)
#axarr[0,1].set_xticks([])

axarr[1,1].plot(age_all_plot, dLi_plot_real, 'g--', label = 'dLi')
#axarr[1].set_ylim(10,20)
#axarr[1,1].set_xlim(0, 5e6)
axarr[1,1].invert_xaxis()
axarr[1,1].legend(loc = 'best')
axarr[1,1].fill_between(age_all_plot, upper_dLi_plot_real, lower_dLi_plot_real, facecolor='green', alpha=0.3)
#axarr[1,1].set_xticks([])

axarr[2,1].plot(age_all_plot, dSr_plot_real, 'b--', label = 'dSr')
axarr[2,1].invert_xaxis()
axarr[2,1].fill_between(age_all_plot, upper_dSr_plot_real, lower_dSr_plot_real, facecolor='blue', alpha=0.3)
#axarr[2].set_ylim(0.7073,0.7078)
axarr[2,1].legend(loc = 'best')


f.tight_layout()
f.subplots_adjust(top = 0.8)

print('number of successes: ' + str(len(success)) + ' & number of failures: ' + str(len(failure)))

#f.subplots_adjust(hspace = 0.5)

#f.savefig('C:/Users/joach/Desktop/Python/Foward_model_figures/4myr_kht1_2.png', bbox_inches='tight')

##Start of filtering
#success = pd.DataFrame()
#failure = pd.DataFrame()
#num_monte_sucess = []
#
#
#if all(dLi_plot <= dLi_plot_real + 0.33*dLi_plot_real) and all(dLi_plot >= dLi_plot_real - 0.33*dLi_plot_real):
#    data = pd.DataFrame({'K_lt': [K_lt_mag], 'K_ht': [K_ht_mag], 'K_riv': [K_riv_mag], 'K_light': [K_light_mag]})
#    success = success.append(data)
#else:
#    data_fail = pd.DataFrame({'K_lt': [K_lt_mag], 'K_ht': [K_ht_mag], 'K_riv': [K_riv_mag], 'K_light': [K_light_mag]})
#    failure = failure.append(data_fail)


#for i in range(0, num_monte):
#    if dLi_plot_real[i] <= dLi_plot_real[i] + 0.33*dLi_plot_real[i] and dLi_plot_real[i] >= dLi_plot_real[i] - 0.33*dLi_plot_real[i]:
#        num_monte_sucess.append(1)
#        if num_monte_success[i] = 1:
#            data = pd.DataFrame({'K_lt': [K_lt_mag], 'K_ht': [K_ht_mag], 'K_riv': [K_riv_mag], 'K_light': [K_light_mag]})
#            success = success.append(data)
#    else:
#        num_monte_success.append(0)
#        data_fail = pd.DataFrame({'K_lt': [-10], 'K_ht': [-10], 'K_riv': [-10], 'K_light': [-10]})
#        success = success.append(data_fail)

#
#
#for i in range(0, len(dLi_plot_real)):
#    if dLi_plot_real[i] <= dLi_plot_real[i] + 0.33*dLi_plot_real[i] and dLi_plot_real[i] >= dLi_plot_real[i] - 0.33*dLi_plot_real[i]:
#        data = pd.DataFrame({'K_lt': [K_lt_mag], 'K_ht': [K_ht_mag], 'K_riv': [K_riv_mag], 'K_light': [K_light_mag]})
#        success = success.append(data)
#    else:
#        data_fail = pd.DataFrame({'K_lt': [-10], 'K_ht': [-10], 'K_riv': [-10], 'K_light': [-10]})
#        success = success.append(data_fail)
#
##
#f, axarr = plt.subplots(3)
#axarr[0].plot(t_plot, dOs_plot, 'r--', label = 'dOs')
##axarr[0].set_ylim(0.4,1)
##axarr[0].plot(test, dOs_ocean_record_smooth[:,1]) this is to plot the real Os record along side the modeled ouput
#axarr[0].legend(loc = 'best')
#axarr[0].set_title(str(len(t_change2)) + ' kyr perturbation w/ ' + 'K_lt: ' + str(K_lt_monte) + ', K_ht: ' + str(K_ht_monte) + ', K_riv: ' + str(K_riv_monte) + ', K_light: ' + str(K_light_monte))
#axarr[0].set_xticks([])
#
#
#axarr[1].plot(t_plot, dLi_plot, 'g--', label = 'dLi')
##axarr[1].set_ylim(10,20)
#axarr[1].legend(loc = 'best')
#axarr[1].set_xticks([])
#
#
#axarr[2].plot(t_plot, dSr_plot, 'b--', label = 'dSr')
##axarr[2].set_ylim(0.7073,0.7078)
#axarr[2].legend(loc = 'best')
#
#f.tight_layout()
#f.subplots_adjust(top = 0.8)
#
#



##dOs_total_origin = dOs_total
###dOs_total = 7.4 * dOs_total/(1 - dOs_total)
##
##dSr_total_origin = dOs_total
##
##dLi_total_origin = dLi_total
##
###Average and std dev of Os outputs
##dOs_ave = np.nanmean(dOs_total, axis = 0)
##dOs_std = np.nanstd(dOs_total, axis = 0)
##dOs_total_low = dOs_ave - dOs_std
##dOs_total_high = dOs_ave + dOs_std
##
###Average and std dev of Sr outputs
##dSr_ave = np.nanmean(dSr_total, axis = 0)
##dSr_std = np.nanstd(dSr_total, axis = 0)
##dSr_total_low = dSr_ave - dSr_std
##dSr_total_high = dSr_ave + dSr_std
##
###Average and std dev of Li outputs
##dLi_ave = np.nanmean(dLi_total, axis = 0)
##dLi_std = np.nanstd(dLi_total, axis = 0)
##dLi_total_low = dLi_ave - dLi_std
##dLi_total_high = dLi_ave + dLi_std
#
#
#age_all_plot_book = np.linspace(AGE_OLD, AGE_YOUNG, 2601)
#
##Swaps columns to rows, and vice versa
#dOs_total_final_t = np.transpose(dOs_total_final)
#dSr_total_final_t = np.transpose(dSr_total_final)
#dLi_total_final_t = np.transpose(dLi_total_final)
#
##Prints all forward modeled outputs to excel
#df_Os = pd.DataFrame(dOs_total_final_t)
#df_Os.insert(0 , 'Age (Ma)', age_all_plot_book)
#df_Os.to_excel('Forwardtest_Os_solvedflux_1Ktoend.xlsx', index = False)
#
#df_Sr = pd.DataFrame(dSr_total_final_t)
#df_Sr.insert(0 , 'Age (Ma)', age_all_plot_book)
#df_Sr.to_excel('Forwardtest_Sr_solvedflux_1Ktoend.xlsx', index = False)
#
#df_Li = pd.DataFrame(dLi_total_final_t)
#df_Li.insert(0 , 'Age (Ma)', age_all_plot_book)
#df_Li.to_excel('Forwardtest_Li_solvedflux_1Ktoend.xlsx', index = False)
#
#
##spreadsheet_data = pd.DataFrame({'Age (Ma)': age_all_plot_book, 'Os_ratio': dOs_total_final_t})
##spreadsheet_data.to_excel('Forwardtest_solvedflux_1Ktoend.xlsx', index = False)
##
#fig, ax = plt.subplots()
#
#ax.plot(t_plot, dOs_plot)
#
#ax.set_ylim(0, 1)

##plt.subplot(2, 1, 1)
##plt.plot(dLi_ocean_record_smooth[:,0], dLi_ocean_record_smooth[:,1])
##plt.ylabel('dLi')
##
##plt.subplot(2,1,2)
##plt.plot(dOs_ocean_record_smooth[:,0], dOs_ocean_record_smooth[:,1])
##
##plt.show()
