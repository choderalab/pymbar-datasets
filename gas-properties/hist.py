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
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.colors import LogNorm
from pylab import *
from matplotlib.ticker import NullFormatter
#=============================================================================================
# PARAMETERS
#=============================================================================================
# Setting temperature array (2 options)
T_option = 2.0 # Set 1.0 or 2.0, according to desired method below

# 1) Automatic, defining initial T and T increment
nominal_Tstart = 300.15 # First NPT input temperature (nominal value) - K
nominal_Tinc = 20.0 # Nominal temperature difference between simulations

# 2) Explicit declared in an array
T_array = [323.15,373.15,473.15,573.15]

# Setting pressure array
nominal_Pstart = [50.0,70.0,90.0] # NPT input pressure (nominal value) - MPa

# Defining states without samples
include_zeros = 0 # Number of points to be included between each isobaric simulation
T_inc = 0.0 # Define temperature diff. between new points. Set "0.0" to equally spaced data based on the first two temperatures
empty_Pstart = [] # Different pressures to be included between each isothermal simulation
#=============================================================================================
# CONSTANTS
#=============================================================================================
kcal2J = 4184.0 # Convert kcal2J
A32m3  = 1E-30 # Convert A3 to m3
Avog = 6.02214129E+23 # Avogadro number
atm2MPa = 0.101325 # Convert atm2MPa
#=============================================================================================
# MAIN
#=============================================================================================
# Print complete arrays (to help debug)
set_printoptions(threshold=nan)

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

temperatureT_H = zeros([1], float64)
pressureT_H = zeros([1], float64)

temperatureT_V = zeros([1], float64)
pressureT_V = zeros([1], float64)

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
        rm_up = n+1+1
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
  range_data = files#(include_zeros+1)*files
  
  pressure = zeros([range_data,T_max], float64)
  temperature = zeros([range_data,T_max], float64)

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
      pressure[ind,k] = x_kt[1,k]*atm2MPa
      temperature[ind,k] = x_kt[3,k]

    if pressureT_H[0] == 0:
      pressureT_H = pressure[ind]
      temperatureT_H = temperature[ind]
      pressureT_V = pressure[ind]
      temperatureT_V = temperature[ind]

    pressureT_H = np.hstack((pressureT_H,pressure[ind]))
    temperatureT_H = np.hstack((temperatureT_H,temperature[ind]))

    pressureT_V = np.vstack((pressureT_V,pressure[ind]))
    temperatureT_V = np.vstack((temperatureT_V,temperature[ind]))
"""
# HEAT MAP
hist2d(temperatureT_H, pressureT_H, bins=50, norm=LogNorm())
colorbar()
show()

# 3D GRAPH
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for k in range((files)*len(temp)+1):

  y = pressureT_V[k]
  x = temperatureT_V[k]

  hist, xedges, yedges = np.histogram2d(x, y, bins=(10,10))

  elements = (len(xedges) - 1) * (len(yedges) - 1)
  xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.25)

  xpos = xpos.flatten()
  ypos = ypos.flatten()
  zpos = np.zeros(elements)
  dx = xedges [1] - xedges [0]
  dy = yedges [1] - yedges [0]
  dz = hist.flatten()

  ax.bar3d(xpos, ypos, zpos, dx, dy, dz, 
             color=np.random.rand(3,1),
             edgecolor='none',
             zsort='average', 
             alpha=1)
ax.set_xlabel('Pressure (MPa)')
ax.set_ylabel('Temperature (K)')
ax.set_zlabel('P(P,T)')
plt.show()
"""
# GRAPH-SET
nbins = 100
filename = 'Histogram'
frameon=False

# Define the x and y data
x = temperatureT_H
y = pressureT_H

# Set up X and Y axes range
xmin = min(x)
xmax = max(x)
ymin = min(y)
ymax = max(y)

# Define the locations for the axes
left, width = 0.12, 0.55
bottom, height = 0.12, 0.55
bottom_h = left_h = left+width+0.02

# Set up the geometry of the three plots
rect_temperature = [left, bottom, width, height] # dimensions of temp plot
rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram

# Set up the size of the figure
fig = plt.figure(1, figsize=(9.5,9))

# Make the three plots
axTemperature = plt.axes(rect_temperature) # temperature plot
axHistx = plt.axes(rect_histx) # x histogram
axHisty = plt.axes(rect_histy) # y histogram

# Remove the inner axes numbers of the histograms
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)
axHistx.yaxis.set_major_formatter(nullfmt)
axHisty.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# Make the 'main' temperature plot
H, xedges,yedges = np.histogram2d(y,x,bins=(nbins,nbins))

# Plot the temperature data
cax = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower',aspect='auto', cmap=cm.jet)
	
#Set up the plot limits
axTemperature.set_xlim(xmin,xmax)
axTemperature.set_ylim(ymin,ymax)

#Set up the histogram limits
axHistx.set_xlim(xmin,xmax)
axHisty.set_ylim(ymin,ymax)

#Plot the histograms
for k in range((files)*len(temp)+1):

  y = pressureT_V[k]
  x = temperatureT_V[k]
  cl = np.random.rand(3,1)

  axHistx.hist(x, bins=nbins/5, color=cl,normed=True,alpha=0.9)
  axHisty.hist(y, bins=nbins/5, orientation='horizontal',color=cl,normed=True,alpha=0.9)

#Plot the axes labels
axTemperature.set_xlabel(r'$\mathit{T}$'+' / K',fontsize=18)
axTemperature.set_ylabel(r'$\mathit{p}$'+ ' / MPa',fontsize=18)
axHistx.set_ylabel('P('+r'$\mathit{T}$'+')',fontsize=18)
axHisty.set_xlabel('P('+r'$\mathit{p}$'+')',fontsize=18, labelpad=20)
axHisty.xaxis.set_label_position("top")

if(filename):
	savefig(filename + '.eps',format = 'eps', transparent=True)
plt.show()

print "DONE!"
