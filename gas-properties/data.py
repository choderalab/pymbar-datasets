# aimoli@gmail.com
# This script extracts data from LAMMPS log files and build the input files to be used as inputs of the run_mbar.py script.
#=============================================================================================
# IMPORTS
#=============================================================================================
import numpy
from numpy import *
from math import *
import pymbar
import timeseries
import commands
import os
import os.path
import sys
import time
import shutil
#=============================================================================================
# PARAMETERS
#=============================================================================================
# Setting temperature array (2 options)
T_option = 1.0 # Set 1.0 or 2.0, according to desired method below

# 1) Automatic, defining initial T and T increment
nominal_Tstart = 300.15 # First NPT input temperature (nominal value) - K
nominal_Tinc = 20.0 # Nominal temperature difference between simulations

# 2) Explicit declared in an array
T_array = []

# Setting pressure array
nominal_Pstart = [50.0,70.0,90.0] # NPT input pressure (nominal value) - MPa

# Defining states without samples
include_zeros = 1 # Number of points to be included between each isobaric simulation
T_inc = 0.0 # Define temperature diff. between new points. Set "0.0" to equally spaced data based on the first two temperatures
empty_Pstart = ["60.0","80.0"] # Different pressures to be included between each isothermal simulation
#=============================================================================================
# CONSTANTS
#=============================================================================================
kcal2J = 4184.0 # Convert kcal2J
A32m3  = 1E-30 # Convert A3 to m3
Avog = 6.02214129E+23 # Avogadro number
#=============================================================================================
# MAIN
#=============================================================================================
# Setting up directories based on pressure values
temp = []
processed_data = []
original_data = []
for i in nominal_Pstart:
  t_m = str(i) + '/temp/'
  p_d = str(i) + '/processed_data/'
  o_d = str(i) + '/original_data/'
  temp.append(t_m)
  processed_data.append(p_d)
  original_data.append(o_d)  

# Start loop over all pressures
Pind = 0
for input_directory in temp:

  if os.path.exists(input_directory):
    shutil.rmtree(input_directory)
    os.mkdir(input_directory)
  else:
    os.mkdir(input_directory)

  nominal_P = nominal_Pstart[Pind]
  output_directory = processed_data[Pind]
  original = original_data[Pind]
  Pind = Pind + 1

  # Defining the limits to extract LAMMPS file data
  for n,line in enumerate(open(original + 'LAMMPS.1','r')):
      if "Step Press" in line: 
        rm_up = n+1
  for m,line in enumerate(open(original + 'LAMMPS.1', 'r')):
      if "Loop time" in line:
        rm_dw = m

  # Extract data from original LAMMPS output files
  files = len(os.walk(original).next()[2]) # Number of files with results to be used
  for ind in range(files):
    # Collect data from files
    filename = os.path.join(original, 'LAMMPS.' + str(ind+1))
    lines = open(filename).readlines()
    open(input_directory + str(ind) + '.lammps', 'w').writelines(lines[rm_up:rm_dw])

  # Define arrays
  filename = os.path.join(input_directory, '0.lammps')
  T_max = int(commands.getoutput('wc -l %s' % filename).split()[0])
  files = len(os.walk(input_directory).next()[2]) # Number of files with results to be used
  range_data = (include_zeros+1)*files
  
  pressure = zeros([range_data], float64)
  pressure_w0s = zeros([range_data], float64)

  vol = zeros([range_data,T_max], float64)
  vol_w0s = zeros([range_data,T_max], float64)

  temperature = zeros([range_data], float64)
  temperature_w0s = zeros([range_data], float64)

  Hconf = zeros([range_data,T_max], float64)
  Hconf_w0s = zeros([range_data,T_max], float64)

  Uconf = zeros([range_data,T_max], float64)
  uconf_w0s = zeros([range_data,T_max], float64)

  print "Reading LAMMPS files"
  for ind in range(files):
    # Collect data from files
    filename = os.path.join(input_directory, str(ind) + '.lammps')
    infile = open(filename, 'r')
    elements = infile.readline().split()
    K = len(elements)

    # Determine maximum number of snapshots
    filename = os.path.join(input_directory, str(ind) + '.lammps')

    # Allocate storage for original energies
    T_k = zeros([K], int32) # T_k[k] is the number of snapshots from umbrella simulation k
    x_kt = zeros([K,T_max], float64) # x_kt[k,t] is the position of snapshot t from energy k (in kcal/mol)

    # Read the energies.
    filename = os.path.join(input_directory, str(ind) + '.lammps')
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    # Parse data.
    for line in lines:
        elements = line.split()
        for k in range(K):
            t = T_k[k]
            x_kt[k,t] = float(elements[k])
            T_k[k] += 1

    for k in range(T_max):
      vol[ind,k] = x_kt[2,k]
      Uconf[ind,k] = x_kt[4,k] - x_kt[16,k] # Calculates the Configurational Internal Energy (PotEn - Emol)
      Hconf[ind,k] = x_kt[4,k] - x_kt[16,k] + nominal_P*1E6*x_kt[2,k]*A32m3/kcal2J*Avog # Calculates the Configurational Enthalpy (PotEn - Emol + Pnominal*Vinstantaneous)

    pressure[ind] = nominal_P
    if T_option == 1.0:
      temperature[ind] = nominal_Tstart + ind * nominal_Tinc
    if T_option == 2.0:
      temperature[ind] = T_array[ind]

  # Generate temperature array
  for i in range(range_data):
    pressure_w0s[i] = nominal_P
  for ind in range(files):
    temperature_w0s[(include_zeros+1)*ind] = temperature[ind]
    for k in range(T_max):
      vol_w0s[(include_zeros+1)*ind,k] = vol[ind,k]
      Hconf_w0s[(include_zeros+1)*ind,k] = Hconf[ind,k]
      uconf_w0s[(include_zeros+1)*ind,k] = Uconf[ind,k]

  if T_inc == 0:
     T_inc = (temperature[1] - temperature[0])/(include_zeros+1)
   
  for ind in range(range_data):
    if temperature_w0s[ind] == 0:
      temperature_w0s[ind] = temperature_w0s[ind-1]+T_inc

  # Write pressures.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'pressure.dat')
     print "Writing Pressure to '%s'..." % filename
     outfile = open(filename, 'w')
     for k in range(0,range_data):
       outfile.write('%.3f' % pressure_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'pressure.dat')
     print "Writing Pressure to '%s'..." % filename
     outfile = open(filename, 'w')
     for k in range(0,range_data-include_zeros):
       outfile.write('%.3f' % pressure_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()

  # Write volumes.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'volume.dat')
     print "Writing Volume to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.4f' % vol_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'volume.dat')
     print "Writing Volume to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.4f' % vol_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()

  # Write temperatures.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'temperature.dat')
     print "Writing Temperature to '%s'..." % filename
     outfile = open(filename, 'w')
     for k in range(0,range_data):
       outfile.write('%.3f' % temperature_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'temperature.dat')
     print "Writing Temperature to '%s'..." % filename
     outfile = open(filename, 'w')
     for k in range(0,range_data-include_zeros):
       outfile.write('%.3f' % temperature_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()

  # Write Hconf.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'hconf.dat')
     print "Writing Hconf to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.6f' % Hconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close() 
  else:
     filename = os.path.join(output_directory, 'hconf.dat')
     print "Writing Hconf to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.6f' % Hconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()  

  # Write Uconf.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'uconf.dat')
     print "Writing Uconf to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.6f' % uconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'uconf.dat')
     print "Writing Uconf to '%s'..." % filename
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.6f' % uconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
  shutil.rmtree(input_directory)

# Create data for new pressures without samples
print "Generating tables of zeros for unsampled pressures (if any)"
processed_data = []
for i in empty_Pstart:
  p_d = str(i) + '/processed_data/'
  processed_data.append(p_d)

# Start loop over all new pressures
Pind=0
for input_directory in empty_Pstart:
  if os.path.exists(input_directory):
    shutil.rmtree(input_directory + '/processed_data/')
    shutil.rmtree(input_directory)
    os.mkdir(input_directory)
    os.mkdir(input_directory + '/processed_data/')
  else:
    os.mkdir(input_directory)
    os.mkdir(input_directory + '/processed_data/')

  nominal_P = empty_Pstart[Pind]
  output_directory = processed_data[Pind]
  Pind = Pind + 1

  pressure = zeros([range_data], float64)
  pressure_w0s = zeros([range_data], float64)

  vol = zeros([range_data,T_max], float64)
  vol_w0s = zeros([range_data,T_max], float64)

  temperature = zeros([range_data], float64)
  temperature_w0s = zeros([range_data], float64)

  Hconf = zeros([range_data,T_max], float64)
  Hconf_w0s = zeros([range_data,T_max], float64)

  Uconf = zeros([range_data,T_max], float64)
  uconf_w0s = zeros([range_data,T_max], float64)

  for ind in range(files):
    pressure[ind] = nominal_P
    temperature[ind] = nominal_Tstart + ind * nominal_Tinc

  # Generate temperature array
  for i in range(range_data):
    pressure_w0s[i] = nominal_P
  for ind in range(files):
    temperature_w0s[(include_zeros+1)*ind] = temperature[ind]

  if T_inc == 0:
     T_inc = (temperature[1] - temperature[0])/(include_zeros+1)
   
  for ind in range(range_data):
    if temperature_w0s[ind] == 0:
      temperature_w0s[ind] = temperature_w0s[ind-1]+T_inc

  # Write new pressures.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'pressure.dat')
     outfile = open(filename, 'w')
     for k in range(0,range_data):
       outfile.write('%.3f' % pressure_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'pressure.dat')
     outfile = open(filename, 'w')
     for k in range(0,range_data-include_zeros):
       outfile.write('%.3f' % pressure_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()

  # Write new volumes.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'volume.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.0f' % vol_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'volume.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.0f' % vol_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()

  # Write new temperatures.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'temperature.dat')
     outfile = open(filename, 'w')
     for k in range(0,range_data):
       outfile.write('%.3f' % temperature_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'temperature.dat')
     outfile = open(filename, 'w')
     for k in range(0,range_data-include_zeros):
       outfile.write('%.3f' % temperature_w0s[k])
       outfile.write(' ')
     outfile.write('\n')
     outfile.close()

  # Write new Hconf.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'hconf.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.0f' % Hconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close() 
  else:
     filename = os.path.join(output_directory, 'hconf.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.0f' % Hconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()  

  # Write new Uconf.
  if include_zeros == 0:
     filename = os.path.join(output_directory, 'uconf.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data):
         outfile.write('%.0f' % uconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
  else:
     filename = os.path.join(output_directory, 'uconf.dat')
     outfile = open(filename, 'w')
     for t in range(T_max):  
       for k in range(0,range_data-include_zeros):
         outfile.write('%.0f' % uconf_w0s[k,t])
         outfile.write(' ')
       outfile.write('\n')
     outfile.close()
print "DONE!"
