# aimoli@gmail.com
# This script uses pyMBAR package to to calculate MBAR estimation of thermodynamic properties via ensemble fluctuations in NPT ensemble.
# The standard estimation is also calculated and the uncertainties are calculated using error propagation equations.
# The ideal contribution to the heat capacity is calculated using the experimental correlation described in U. Setzmann, W. Wagner, J. Phys. Chem. Ref. Data, 20 (1991) 1061-1150
#=============================================================================================
# IMPORTS
#=============================================================================================
# Requires the Uncertainties package available at https://pypi.python.org/pypi/uncertainties/
from uncertainties import unumpy
import numpy
import numpy as np
from numpy import *
from math import *
import pymbar
import timeseries
import commands
import os
import os.path
import sys
import time
#=============================================================================================
# PARAMETERS
#=============================================================================================
correlated_data = 1  # Set "0" for uncorrelated original data and "1" to use subsampling
nominal_Pstart = [50.0,70.0,90.0] # NPT input pressure (nominal value) - MPa
N_CH4 = 512 # Number of CH4 molecules on the simulations
N_CO2 = 0 # Number of CO2 molecules on the simulations
N_total = N_CO2 + N_CH4 # Total Number of molecules on the simulations
#=============================================================================================
# CONSTANTS
#=============================================================================================
kcal2J = 4184.0 # Convert kcal2J
A32m3  = 1E-30 # Convert A3 to m3
Cmass  = 0.0120107 # Carbon mass (kg/mol)
Omass = 0.0159994 # Oxigen mass (kg/mol)
Hmass = 0.00100794 # Hidrogen mass (kg/mol)
Avog = 6.02214129E+23 # Avogadro number
kB_JK = 1.3806488E-23 # Boltzmann constant in J/K
kB_kcalmolK = kB_JK * Avog / kcal2J # Boltzmann constant in kcal/mol/K
Tc = 304.1282 # Span and Wagner Cp correlation parameteres
a0 = [0.0,8.37304456,-3.70454304,2.50000000,1.99427042,0.62105248,0.41195293,1.04028922,0.08327678] # Span and Wagner Cp correlation parameteres
theta0 = [0.0,0.0,0.0,0.0,3.15163,6.11190,6.77708,11.32384,27.08792] # Span and Wagner Cp correlation parameteres
n = [4.0016,0.008449,4.6942,3.4865,1.6572,1.4115] # Setzmann Cp correlation parameteres
theta = [0.0,648,1957,3895,5705,15080] # Setzmann Cp correlation parameteres
#=============================================================================================
# MAIN
#=============================================================================================
# Print complete arrays (to help debug)
set_printoptions(threshold=nan)

# Setting up directories based on pressure values
processed_data = []
mbar_results = []
plots = []
NIST_data = []
for i in nominal_Pstart:
    p_d = str(i) + '/processed_data/'
    m_r = str(i) + '/mbar_results/'
    pl_d = m_r + 'plots/'
    n_d = str(i) + '/NIST_data/'
    processed_data.append(p_d)
    mbar_results.append(m_r)
    plots.append(pl_d)
    NIST_data.append(n_d)

# Start loop over all pressure
Pind = 0
for directory in processed_data:
    output_directory = mbar_results[Pind]
    plot_directory = plots[Pind]
    NIST_directory = NIST_data[Pind]
    Pind = Pind + 1
    
    print "Reading input files from %s directory..." % directory
    # Read Temperature
    filename = os.path.join(directory,'temperature.dat')
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    elements = infile.readline().split()
    K = len(elements)
    temperature = zeros([K], float64)
    for k in range(K):
        temperature[k] = float(elements[k])
    infile.close()
    
    # Read Pressure
    filename = os.path.join(directory,'pressure.dat')
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    elements = infile.readline().split()
    pressure = zeros([K], float64)
    for k in range(K):
        pressure[k] = float(elements[k])
    infile.close()
    
    # Determine maximum number of snapshots in all simulations
    filename = os.path.join(directory,'hconf.dat')
    T_max = int(commands.getoutput('wc -l %s' % filename).split()[0])
    
    # Allocate storage for original Hconf
    T_k = zeros([K], int32)
    hconf_original = zeros([K,T_max], float64)
    
    # Read Hconf
    filename = os.path.join(directory,'hconf.dat')
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    
    # Parse data
    for line in lines:
        elements = line.split()
        for k in range(K):
            t = T_k[k]
            hconf_original[k,t] = float(elements[k])
            T_k[k] += 1
    
    # Allocate storage for original Volume
    T_k = zeros([K], int32)
    volume_original = zeros([K,T_max], float64)
    
    # Read Volume
    filename = os.path.join(directory,'volume.dat')
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    
    # Parse data
    for line in lines:
        elements = line.split()
        for k in range(K):
            t = T_k[k]
            volume_original[k,t] = float(elements[k])
            T_k[k] += 1

    # Allocate storage for original Uconf
    T_k = zeros([K], int32)
    uconf_original = zeros([K,T_max], float64)

    # Read Uconf
    filename = os.path.join(directory,'uconf.dat')
    print "Reading %s..." % filename
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    # Parse data
    for line in lines:
        elements = line.split()
        for k in range(K):
            t = T_k[k]
            uconf_original[k,t] = float(elements[k])
            T_k[k] += 1

    N_max = 0
    for k in range(K):
        N_max = max(N_max, ceil(T_k[k]))

    # Subsample trajectories based on Hconf
    hconf = zeros([K, N_max], float64)
    volume = zeros([K, N_max], float64)
    uconf = zeros([K, N_max], float64)
    N_k = zeros([K], int32)

    if correlated_data == 1:
        hconf = zeros([K, T_max], float64)
        N_ksam = zeros([K], int32)
        indices2 = zeros([T_max], int32)
        for k in range(1,T_max):
            indices2[k] = k
        for k in range(K):
            # Compute correlation times
            indices = timeseries.subsampleCorrelatedData(hconf_original[k,0:N_k[k]])
            # Store subsampled positions
            if len(indices) >= 1000:
                N_ksam[k] = len(indices)
                hconf[k,0:N_ksam[k]] = hconf_original[k,indices]
                volume[k,0:N_ksam[k]] = volume_original[k,indices]
                uconf[k,0:N_ksam[k]] = uconf_original[k,indices]
            else:
                N_ksam[k] = len(indices2)
                hconf[k,0:N_ksam[k]] = hconf_original[k,indices2]
                volume[k,0:N_ksam[k]] = volume_original[k,indices2]
                uconf[k,0:N_ksam[k]] = uconf_original[k,indices2]
            print('\n')
        N_k = zeros([K], int32)
        N_k = N_ksam

    for k in range(len(temperature)):
        if hconf[k,0] == 0:
            N_k[k] = 0
        else:
            N_k[k] = T_max
    #=============================================================================================
    # RUNNING MBAR
    #=============================================================================================
    # Calculate reduced potentials
    beta_k = (kB_kcalmolK * temperature)**(-1)
    u_kln = zeros([K,K,N_max], float64)
    hconfm = zeros([K,K,N_max], float64)
    vol_hconfm = zeros([K,K,N_max], float64)
    uconf_hconfm = zeros([K,K,N_max], float64)

    for k in range(K):
        for l in range(K):
            hconfm[k,l,0:N_k[k]] = uconf[k,0:N_k[k]] + pressure[l]*1E6*volume[k,0:N_k[k]]*A32m3/kcal2J*Avog
            u_kln[k,l,0:N_k[k]] = beta_k[l] * hconfm[k,l,0:N_k[k]]
            vol_hconfm[k,l,0:N_k[k]] = volume[k,0:N_k[k]] * hconfm[k,l,0:N_k[k]]
            uconf_hconfm[k,l,0:N_k[k]] = uconf[k,0:N_k[k]] * hconfm[k,l,0:N_k[k]]

    # Read initial guess for free energies (if there is a file)
    if os.path.exists(output_directory + 'f_k.dat'):
        print " "
        print "Reading free energies from f_k.dat"
        infile = open(output_directory + 'f_k.dat', 'r')
        lines = infile.readlines()
        infile.close()
        elements = list()
        for line in lines:
            elements += line.split()
        K = len(elements)
        f_k = numpy.zeros([K], numpy.float64)
        for k in range(K):
            f_k[k] = float(elements[k])
        
        # Initialize MBAR
        print "Running MBAR..."
        mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', relative_tolerance = 1.0e-10, initial_f_k = f_k)
    else:
        print "Running MBAR..."
        mbar = pymbar.MBAR(u_kln, N_k, verbose = True, method = 'adaptive', relative_tolerance = 1.0e-10)
    #=============================================================================================
    # COMPUTING OBSERVABLES
    #=============================================================================================
    # Preparing arrays for standard observables calculation (removing columns of zeros, if any)
    volume2 = numpy.zeros([K,N_max], numpy.float32)
    vol_hconf = numpy.zeros([K,N_max], numpy.float32)
    uconf_hconf = numpy.zeros([K,N_max], numpy.float32)

    for k in range(K):
       volume2[k,0:N_k[k]] = volume[k,0:N_k[k]]*volume[k,0:N_k[k]]
       vol_hconf[k,0:N_k[k]] = volume[k,0:N_k[k]]*hconf[k,0:N_k[k]]
       uconf_hconf[k,0:N_k[k]] = uconf[k,0:N_k[k]]*hconf[k,0:N_k[k]]

    i = 0
    for k in range(K):
        if not hconf[k,0] == 0:
            i += 1

    N_k_std = numpy.zeros([i], numpy.float64)
    temperature_std = numpy.zeros([i], numpy.float64)
    pressure_std = numpy.zeros([i], numpy.float64)
    hconf_std = numpy.zeros([i, N_max], numpy.float64)
    volume_std = numpy.zeros([i, N_max], numpy.float64)
    volume2_std = numpy.zeros([i, N_max], numpy.float64)
    vol_hconf_std = numpy.zeros([i, N_max], numpy.float64)
    uconf_std = numpy.zeros([i, N_max], numpy.float64)
    uconf_hconf_std = numpy.zeros([i, N_max], numpy.float64)

    hconf_STD = numpy.zeros([i],numpy.float64)
    dhconf_STD = numpy.zeros([i],numpy.float64)

    volume_STD = numpy.zeros([i],numpy.float64)
    dvolume_STD = numpy.zeros([i],numpy.float64)

    volume2_STD = numpy.zeros([i],numpy.float64)
    dvolume2_STD = numpy.zeros([i],numpy.float64)

    vol_hconf_STD = numpy.zeros([i],numpy.float64)
    dvol_hconf_STD = numpy.zeros([i],numpy.float64)

    uconf_STD = numpy.zeros([i],numpy.float64)
    duconf_STD = numpy.zeros([i],numpy.float64)

    uconf_hconf_STD = numpy.zeros([i],numpy.float64)
    duconf_hconf_STD = numpy.zeros([i],numpy.float64)

    # Removing columns of zeros, if any
    I = 0
    for k in range(K):
        if not hconf[k,0] == 0:
            hconf_std[I] = hconf[k]
            N_k_std[I] = N_k[k]
            temperature_std[I] = temperature[k]
            pressure_std[I] = pressure[k]
            volume_std[I] = volume[k]
            volume2_std[I] = volume2[k]
            vol_hconf_std[I] = vol_hconf[k]
            uconf_std[I] = uconf[k]
            uconf_hconf_std[I] = uconf_hconf[k]
            I += 1

    # Calculating standard estimates of observables and deviations
    for k in range(I):
        hconf_STD[k] = numpy.average(hconf_std[k,0:N_k_std[k]])
        dhconf_STD[k]  = numpy.sqrt(numpy.var(hconf_std[k,0:N_k_std[k]])/(N_k_std[k]-1))
        volume_STD[k] = numpy.average(volume_std[k,0:N_k_std[k]])
        dvolume_STD[k]  = numpy.sqrt(numpy.var(volume_std[k,0:N_k_std[k]])/(N_k_std[k]-1))
        volume2_STD[k] = numpy.average(volume2_std[k,0:N_k_std[k]])
        dvolume2_STD[k]  = numpy.sqrt(numpy.var(volume2_std[k,0:N_k_std[k]])/(N_k_std[k]-1))
        vol_hconf_STD[k] = numpy.average(vol_hconf_std[k,0:N_k_std[k]])
        dvol_hconf_STD[k]  = numpy.sqrt(numpy.var(vol_hconf_std[k,0:N_k_std[k]])/(N_k_std[k]-1))
        uconf_STD[k] = numpy.average(uconf_std[k,0:N_k_std[k]])
        duconf_STD[k]  = numpy.sqrt(numpy.var(uconf_std[k,0:N_k_std[k]])/(N_k_std[k]-1))
        uconf_hconf_STD[k] = numpy.average(uconf_hconf_std[k,0:N_k_std[k]])
        duconf_hconf_STD[k]  = numpy.sqrt(numpy.var(uconf_hconf_std[k,0:N_k_std[k]])/(N_k_std[k]-1))

    # Using MBAR to estimate observables
    print "Calculating observables"
    print "Hconf..."
    (hconf_MBAR, dhconf_MBAR) = mbar.computeExpectations(hconfm)
    print "Volume..."
    (volume_MBAR, dvolume_MBAR) = mbar.computeExpectations(volume)
    print "Volume^2..."
    (volume2_MBAR, dvolume2_MBAR) = mbar.computeExpectations(volume2)
    print "Volume x Hconf..."
    (vol_hconf_MBAR, dvol_hconf_MBAR) = mbar.computeExpectations(vol_hconfm)
    print "Uconf..."
    (uconf_MBAR, duconf_MBAR) = mbar.computeExpectations(uconf)
    print "Uconf x Hconf..."
    (uconf_hconf_MBAR, duconf_hconf_MBAR) = mbar.computeExpectations(uconf_hconfm)
    #=============================================================================================
    # CALCULATING PROPERTIES
    #=============================================================================================
    # Preparing arrays
    temperature_dev = numpy.zeros([K],numpy.float64)
    dtemperature = numpy.zeros([K],numpy.float64)

    pressure_dev = numpy.zeros([K],numpy.float64)
    dpressure = numpy.zeros([K],numpy.float64)

    temperature_std_dev = numpy.zeros([I],numpy.float64)
    dtemperature_std = numpy.zeros([I],numpy.float64)

    pressure_std_dev = numpy.zeros([I],numpy.float64)
    dpressure_std = numpy.zeros([I],numpy.float64)

    hconf_STDdev = numpy.zeros([I],numpy.float64)
    hconf_MBARdev = numpy.zeros([K],numpy.float64)

    volume_STDdev = numpy.zeros([I],numpy.float64)
    volume_MBARdev = numpy.zeros([K],numpy.float64)

    volume2_STDdev = numpy.zeros([I],numpy.float64)
    volume2_MBARdev = numpy.zeros([K],numpy.float64)

    vol_hconf2_STDdev = numpy.zeros([I],numpy.float64)
    vol_hconf2_MBARdev = numpy.zeros([K],numpy.float64)

    uconf_STDdev = numpy.zeros([I],numpy.float64)
    uconf_MBARdev = numpy.zeros([K],numpy.float64)

    uconf_hconf_STDdev = numpy.zeros([I],numpy.float64)
    duconf_hconf_STDdev = numpy.zeros([K],numpy.float64)

    rho_STDdev = numpy.zeros([I],numpy.float64)
    rho_MBARdev = numpy.zeros([K],numpy.float64)

    aP_STDdev = numpy.zeros([I],numpy.float64)
    aP_MBARdev = numpy.zeros([K],numpy.float64)

    kT_STDdev = numpy.zeros([I],numpy.float64)
    kT_MBARdev = numpy.zeros([K],numpy.float64)

    Cp_id_STD = numpy.zeros([I],numpy.float64)
    Cp_id_MBAR = numpy.zeros([K],numpy.float64)

    Cp_STDdev = numpy.zeros([I],numpy.float64)
    Cp_MBARdev = numpy.zeros([K],numpy.float64)

    Cv_STDdev = numpy.zeros([I],numpy.float64)
    Cv_MBARdev = numpy.zeros([K],numpy.float64)

    uJT_STDdev = numpy.zeros([I],numpy.float64)
    uJT_MBARdev = numpy.zeros([K],numpy.float64)

    SS_STDdev = numpy.zeros([I],numpy.float64)
    SS_MBARdev = numpy.zeros([K],numpy.float64)

    # Concatenate averages and deviations into the same array
    temperature_dev = unumpy.uarray((temperature,dtemperature))

    temperature_std_dev = unumpy.uarray((temperature_std,dtemperature_std))

    pressure_dev = unumpy.uarray((pressure,dpressure))

    pressure_std_dev = unumpy.uarray((pressure_std,dpressure_std))

    hconf_STDdev = unumpy.uarray((hconf_STD,dhconf_STD))
    hconf_MBARdev = unumpy.uarray((hconf_MBAR,dhconf_MBAR))

    hconf_STDdev = unumpy.uarray((hconf_STD,dhconf_STD))
    hconf_MBARdev = unumpy.uarray((hconf_MBAR,dhconf_MBAR))

    volume_STDdev = unumpy.uarray((volume_STD,dvolume_STD))
    volume_MBARdev = unumpy.uarray((volume_MBAR,dvolume_MBAR))

    volume2_STDdev = unumpy.uarray((volume2_STD,dvolume2_STD))
    volume2_MBARdev = unumpy.uarray((volume2_MBAR,dvolume2_MBAR))

    vol_hconf_STDdev = unumpy.uarray((vol_hconf_STD,dvol_hconf_STD))
    vol_hconf_MBARdev = unumpy.uarray((vol_hconf_MBAR,dvol_hconf_MBAR))

    uconf_STDdev = unumpy.uarray((uconf_STD,duconf_STD))
    uconf_MBARdev = unumpy.uarray((uconf_MBAR,duconf_MBAR))

    uconf_hconf_STDdev = unumpy.uarray((uconf_hconf_STD,duconf_hconf_STD))
    uconf_hconf_MBARdev = unumpy.uarray((uconf_hconf_MBAR,duconf_hconf_MBAR))

    # rho (kg/m3)
    rho_STDdev = ((Cmass+2*Omass)*N_CO2+(Cmass+4*Hmass)*N_CH4)/Avog/(volume_STDdev*1E-30)
    rho_MBARdev = ((Cmass+2*Omass)*N_CO2+(Cmass+4*Hmass)*N_CH4)/Avog/(volume_MBARdev*1E-30)

    # aP (1/K)
    aP_STDdev = (vol_hconf_STDdev-volume_STDdev*hconf_STDdev)/(volume_STDdev*kB_kcalmolK*temperature_std_dev*temperature_std_dev)
    aP_MBARdev = (vol_hconf_MBARdev-volume_MBARdev*hconf_MBARdev)/(volume_MBARdev*kB_kcalmolK*temperature_dev*temperature_dev)

    # kT (1/Pa)
    kT_STDdev = (volume2_STDdev-volume_STDdev*volume_STDdev)/(volume_STDdev*kB_JK*temperature_std_dev)*A32m3
    kT_MBARdev = (volume2_MBARdev-volume_MBARdev*volume_MBARdev)/(volume_MBARdev*kB_JK*temperature_dev)*A32m3

    # Cp_id - Experimental correlations from Setzmann for CH4 (J/K)
    for k in range(I):
        Cp_id_STD[k] = N_CO2*(kB_JK*(1+a0[3]+a0[4]*(theta0[4]*Tc/temperature_std[k])**2*(exp(theta0[4]*Tc/temperature_std[k]))/(exp(theta0[4]*Tc/temperature_std[k])-1)**2+a0[5]*(theta0[5]*Tc/temperature_std[k])**2*(exp(theta0[5]*Tc/temperature_std[k]))/(exp(theta0[5]*Tc/temperature_std[k])-1)**2+a0[6]*(theta0[6]*Tc/temperature_std[k])**2*(exp(theta0[6]*Tc/temperature_std[k]))/(exp(theta0[6]*Tc/temperature_std[k])-1)**2+a0[7]*(theta0[7]*Tc/temperature_std[k])**2*(exp(theta0[7]*Tc/temperature_std[k]))/(exp(theta0[7]*Tc/temperature_std[k])-1)**2+a0[8]*(theta0[8]*Tc/temperature_std[k])**2*(exp(theta0[8]*Tc/temperature_std[k]))/(exp(theta0[8]*Tc/temperature_std[k])-1)**2))+ N_CH4*(kB_JK*(n[0]+n[1]*(theta[1]/temperature_std[k])**2*exp(theta[1]/temperature_std[k])/(exp(theta[1]/temperature_std[k])-1)**2+n[2]*(theta[2]/temperature_std[k])**2*exp(theta[2]/temperature_std[k])/(exp(theta[2]/temperature_std[k])-1)**2+n[3]*(theta[3]/temperature_std[k])**2*exp(theta[3]/temperature_std[k])/(exp(theta[3]/temperature_std[k])-1)**2+n[4]*(theta[4]/temperature_std[k])**2*exp(theta[4]/temperature_std[k])/(exp(theta[4]/temperature_std[k])-1)**2+n[5]*(theta[5]/temperature_std[k])**2*exp(theta[5]/temperature_std[k])/(exp(theta[5]/temperature_std[k])-1)**2))

    for k in range(K):
        Cp_id_MBAR[k] = N_CO2*(kB_JK*(1+a0[3]+a0[4]*(theta0[4]*Tc/temperature[k])**2*(exp(theta0[4]*Tc/temperature[k]))/(exp(theta0[4]*Tc/temperature[k])-1)**2+a0[5]*(theta0[5]*Tc/temperature[k])**2*(exp(theta0[5]*Tc/temperature[k]))/(exp(theta0[5]*Tc/temperature[k])-1)**2+a0[6]*(theta0[6]*Tc/temperature[k])**2*(exp(theta0[6]*Tc/temperature[k]))/(exp(theta0[6]*Tc/temperature[k])-1)**2+a0[7]*(theta0[7]*Tc/temperature[k])**2*(exp(theta0[7]*Tc/temperature[k]))/(exp(theta0[7]*Tc/temperature[k])-1)**2+a0[8]*(theta0[8]*Tc/temperature[k])**2*(exp(theta0[8]*Tc/temperature[k]))/(exp(theta0[8]*Tc/temperature[k])-1)**2))+ N_CH4*(kB_JK*(n[0]+n[1]*(theta[1]/temperature[k])**2*exp(theta[1]/temperature[k])/(exp(theta[1]/temperature[k])-1)**2+n[2]*(theta[2]/temperature[k])**2*exp(theta[2]/temperature[k])/(exp(theta[2]/temperature[k])-1)**2+n[3]*(theta[3]/temperature[k])**2*exp(theta[3]/temperature[k])/(exp(theta[3]/temperature[k])-1)**2+n[4]*(theta[4]/temperature[k])**2*exp(theta[4]/temperature[k])/(exp(theta[4]/temperature[k])-1)**2+n[5]*(theta[5]/temperature[k])**2*exp(theta[5]/temperature[k])/(exp(theta[5]/temperature[k])-1)**2))

    # Cp (J/K)
    Cp_STDdev = ((uconf_hconf_STDdev-uconf_STDdev*hconf_STDdev)*(kcal2J/Avog)**2/(kB_JK*temperature_std_dev*temperature_std_dev))+((vol_hconf_STDdev-volume_STDdev*hconf_STDdev)*(kcal2J/Avog*A32m3)*(pressure_std_dev*1E6)/(kB_JK*temperature_std_dev*temperature_std_dev))-(N_total*kB_JK) + Cp_id_STD
    Cp_MBARdev = ((uconf_hconf_MBARdev-uconf_MBARdev*hconf_MBARdev)*(kcal2J/Avog)**2/(kB_JK*temperature_dev*temperature_dev))+((vol_hconf_MBARdev-volume_MBARdev*hconf_MBARdev)*(kcal2J/Avog*A32m3)*(pressure_dev*1E6)/(kB_JK*temperature_dev*temperature_dev))-(N_total*kB_JK) + Cp_id_MBAR

    # Cv (1/Pa)
    Cv_STDdev = (Cp_STDdev-temperature_std_dev*(volume_STDdev*A32m3)*aP_STDdev**2/kT_STDdev)
    Cv_MBARdev = (Cp_MBARdev-temperature_dev*(volume_MBARdev*A32m3)*aP_MBARdev**2/kT_MBARdev)

    # uJT (K/Pa)
    uJT_STDdev = (volume_STDdev*A32m3)/Cp_STDdev*(temperature_std_dev*aP_STDdev-1)
    uJT_MBARdev = (volume_MBARdev*A32m3)/Cp_MBARdev*(temperature_dev*aP_MBARdev-1)

    # SS (m/s) - Sets the value to zero in case of unphysical behaviour that might lead to sqrt of negative numbers
    SS_STDdev = (Cp_STDdev/Cv_STDdev/kT_STDdev/rho_STDdev)
    SS_MBARdev = (Cp_MBARdev/Cv_MBARdev/kT_MBARdev/rho_MBARdev)

    for k in range(I):
        if unumpy.nominal_values(SS_STDdev[k]) >= 0.0:
            SS_STDdev[k] = SS_STDdev[k]**.5
        else:
            SS_STDdev[k] = 0.0

    for k in range(K):
        if unumpy.nominal_values(SS_MBARdev[k]) >= 0.0:
            SS_MBARdev[k] = SS_MBARdev[k]**.5
        else:
            SS_MBARdev[k] = 0.0

    # Convert units
    kT_STDdev = kT_STDdev*1E6
    kT_MBARdev = kT_MBARdev*1E6

    Cp_id_STD = Cp_id_STD*Avog/N_total
    Cp_id_MBAR = Cp_id_MBAR*Avog/N_total

    Cp_STDdev = Cp_STDdev*Avog/N_total
    Cp_MBARdev = Cp_MBARdev*Avog/N_total

    Cv_STDdev = Cv_STDdev*Avog/N_total
    Cv_MBARdev = Cv_MBARdev*Avog/N_total

    uJT_STDdev = uJT_STDdev*1E6
    uJT_MBARdev = uJT_MBARdev*1E6

    #=============================================================================================
    # OUTPUT RESULTS
    #=============================================================================================
    # Write observables
    filename = os.path.join(output_directory,'STD_results.dat')
    print "Writing final results to '%s'..." % filename
    outfile = open(filename, 'w')
    outfile.write('Temperature(K) Pressure(MPa) Hconf(kcal/mol) dHconf Volume(A3) dV V^2(A6) dV^2 V*Hconf(A3*kcal/mol) dV*Hconf Uconf(kcal/mol) dUconf Uconf*Hconf(kcal/mol) dUconf*Hconf rho(kg/m3) drho aP(1/K) daP kT(1/MPa) dkT Cp_id(J/molK) Cp(J/molK) dCp Cv(J/molK) dCv uJT(K/MPa) duJT SS(m/s) dSS')
    outfile.write('\n')
    for k in range(I):
        outfile.write('%.10e' % unumpy.nominal_values(temperature_std_dev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(pressure_std_dev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(volume_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(volume_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(volume2_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(volume2_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(vol_hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(vol_hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uconf_hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uconf_hconf_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(rho_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(rho_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(aP_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(aP_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(kT_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(kT_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cp_id_STD[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cp_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(Cp_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cv_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(Cv_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uJT_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uJT_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(SS_STDdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(SS_STDdev[k]))
        outfile.write('\n')
    outfile.close()

    filename = os.path.join(output_directory,'MBAR_results.dat')
    print "Writing final results to '%s'..." % filename
    outfile = open(filename, 'w')
    outfile.write('Temperature(K) Pressure(MPa) Hconf(kcal/mol) dHconf Volume(A3) dV V^2(A6) dV^2 V*Hconf(A3*kcal/mol) dV*Hconf Uconf(kcal/mol) dUconf Uconf*Hconf(kcal/mol) dUconf*Hconf rho(kg/m3) drho aP(1/K) daP kT(1/MPa) dkT Cp_id(J/molK) Cp(J/molK) dCp Cv(J/molK) dCv uJT(K/MPa) duJT SS(m/s) dSS')
    outfile.write('\n')
    for k in range(K):
        outfile.write('%.10e' % unumpy.nominal_values(temperature_dev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(pressure_dev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(volume_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(volume_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(volume2_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(volume2_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(vol_hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(vol_hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uconf_hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uconf_hconf_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(rho_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(rho_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(aP_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(aP_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(kT_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(kT_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cp_id_MBAR[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cp_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(Cp_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(Cv_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(Cv_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(uJT_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(uJT_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.nominal_values(SS_MBARdev[k]))
        outfile.write(' ')
        outfile.write('%.10e' % unumpy.std_devs(SS_MBARdev[k]))
        outfile.write('\n')
    outfile.close()

    # Write free energies
    print "Writing free energies"
    f_k=mbar.f_k
    K = f_k.size
    outfile = open(output_directory + 'f_k.dat', 'w')
    for k in range(K):
        outfile.write("%.10e " % f_k[k])
    outfile.write("\n")
    outfile.close()

    #=============================================================================================
    # GNUPLOT PLOTS
    #=============================================================================================
    # make gnuplot plots
    data_STD = '%(output_directory)sSTD_results.dat' % vars()
    data_MBAR = '%(output_directory)sMBAR_results.dat' % vars()
    NIST1 = '%(NIST_directory)sNIST1.dat' % vars()
    NIST3 = '%(NIST_directory)sNIST3.dat' % vars()
    NIST5 = '%(NIST_directory)sNIST5.dat' % vars()
    gnuplot_in = os.path.join(plot_directory + 'gnuplot.in')

    # Density
    filename = os.path.join(plot_directory, 'rho.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Density"
        set xlabel "Temperature (K)"
        set ylabel "Density (kg/m^3)"
        plot "%(data_STD)s" u 1:15 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:15 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:3 with lines t "NIST", "%(NIST3)s" u 1:3 with lines notitle, "%(NIST5)s" u 1:3 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())

    # Coefficient of Thermal Expansion
    filename = os.path.join(plot_directory, 'aP.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Coefficient of Thermal Expansion"
        set xlabel "Temperature (K)"
        set ylabel "aP (1/K)"
        plot "%(data_STD)s" u 1:17 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:17 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:14 with lines t "NIST", "%(NIST3)s" u 1:14 with lines notitle, "%(NIST5)s" u 1:14 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())

    # Isothermal Compressibility
    filename = os.path.join(plot_directory, 'kT.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Isothermal Compressibility"
        set xlabel "Temperature (K)"
        set ylabel "aP (1/MPa)"
        plot "%(data_STD)s" u 1:19 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:19 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:15 with lines t "NIST", "%(NIST3)s" u 1:15 with lines notitle, "%(NIST5)s" u 1:15 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())

    # Isobaric Heat Capacity
    filename = os.path.join(plot_directory, 'Cp.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Isobaric Heat Capacity"
        set xlabel "Temperature (K)"
        set ylabel "Cp (J/mol.K)"
        plot "%(data_STD)s" u 1:22 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:22 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:9 with lines t "NIST", "%(NIST3)s" u 1:9 with lines notitle, "%(NIST5)s" u 1:9 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())
    
    # Isochoric Heat Capacity
    filename = os.path.join(plot_directory, 'Cv.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Isochoric Heat Capacity"
        set xlabel "Temperature (K)"
        set ylabel "Cv (J/mol.K)"
        plot "%(data_STD)s" u 1:24 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:24 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:8 with lines t "NIST", "%(NIST3)s" u 1:8 with lines notitle, "%(NIST5)s" u 1:8 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())
    
    # Joule-Thomson Coefficient
    filename = os.path.join(plot_directory, 'uJT.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Joule-Thomson Coefficient"
        set xlabel "Temperature (K)"
        set ylabel "uJT (K/MPa)"
        plot "%(data_STD)s" u 1:26 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:26 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:11 with lines t "NIST", "%(NIST3)s" u 1:11 with lines notitle, "%(NIST5)s" u 1:11 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())
    
    # Speed of Sound
    filename = os.path.join(plot_directory, 'SS.eps')
    gnuplot_input = """
        set term postscript solid
        set output "%(filename)s"
        set title "Speed of Sound"
        set xlabel "Temperature (K)"
        set ylabel "SS (m/s)"
        plot "%(data_STD)s" u 1:28 t "STD estimate" with points pt 7 lc rgb "red" ps 1, "%(data_MBAR)s" u 1:28 t "MBAR optimal estimate" with points pt 6 lc rgb "black" ps 1, "%(NIST1)s" u 1:10 with lines t "NIST", "%(NIST3)s" u 1:10 with lines notitle, "%(NIST5)s" u 1:10 with lines notitle
        """ % vars()
    outfile = open(gnuplot_in, 'w')
    outfile.write(gnuplot_input)
    outfile.close()
    output = commands.getoutput('gnuplot < %(gnuplot_in)s' % vars())
    
    output = commands.getoutput('rm %(gnuplot_in)s' % vars())

print "DONE!"
