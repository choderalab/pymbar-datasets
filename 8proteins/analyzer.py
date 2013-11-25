#!/usr/bin/python

#=========================================================
#IMPORTS
#=========================================================
import sys
import numpy
from numpy import *
from math import *
import pymbar # for MBAR analysis
import timeseries # for timeseries analysis
import commands
import os
import os.path
import pdb

#=========================================================
# CONSTANTS
#=========================================================
kB = 0.00831447  #Boltzmann constant (Gas constant) in kJ/(mol*K)
MinT = 200.0  # Minimum of the Temperature Range in K
MaxT = 400.0  # Maximum of the Temperature Range in K
MinTrans = 338  #Start temperature of tranisition (melting point) to allow for 2nd increment  
MaxTrans = 353  #End Temperature of transition
increment1 = 1.5  # Temperature increment for non-transition zone (flat, needs less)
increment2 = 0.4  # Temperature increment for transition zone (dynamic, needs as many points as possible)
#=========================================================
# PARAMETERS
#=========================================================
data_directory = os.path.join('in_mbar')
out_directory = os.path.join('out_mbar')
temperature_list_filename = os.path.join(data_directory, 'temperatures')
potential_energies_filename = os.path.join(data_directory, 'p_e') 
total_energies_filename = os.path.join(data_directory, 'tot_e')
contacts_filename = os.path.join(data_directory, 'contacts')
#-------------Don't Change Things Below Here--------------#

#=========================================================
# SUBROUTINES
#=========================================================

def read_file(filename):
	"""Read contents of the specified file.

	ARGUMENTS
	  filename (string) - the name of the file to be read

	RETURNS
	  lines (list of strings) - the contents of the file, split by line
	"""

	infile = open(filename, 'r')
	lines = infile.readlines()
	infile.close()
	
	return lines

def logSum(log_terms):
	"""Compute the log of a sum whose logarithms are provided.

	ARGUMENTS
	  log_terms is the array (possibly multidimensional) containing the logs of the terms to be summed.
	
	RETURNS
	  log_sum is the log of the sum of the terms
	"""

  #compute the maximum argument
	max_log_term = log_terms.max()
  
  #compute the reduced terms
	terms = numpy.exp(log_terms - max_log_term)

  #compute the log sum
	log_sum = log(terms.sum() ) + max_log_term

  #return the log sum
	return log_sum

def histogram_wham(beta_k, U_kn, N_k, nbins = 100, bind_width = None, maximum_iterations = 5000, relative_tolerance = 1.0e-8, initial_f_k = None):
	"""Construct an initial guess of the f_k by histogram reweighting (specifically, WHAM).
	
	ARGUMENTS
	  beta_k (numpy K array) - inverse temperatures (in units of 1/energy)
	  U_kn (numpy K x N_max array) - potential energies (in units of energy)
	  N_k (numpy K array of int32) - number of samples per states, N_max = N_k.max()

	OPTIONAL ARGUMENTS
	  nbins (int) - if specified, sets the number of bins to use (default = 100)
	  bin_width (float) - if specified, sets the bin width (overrides nbins) (default = None)
	  maximum_iterations (int) - maximum number of iterations to use
	  relative_tolerance (float) - relative convergence tolerance (default = 1.0e-8)

	RETURNS
	  f_k (numpy K array) - guess at initial state dimensionless free energies
	"""

  # Get sizes
	K = N_k.size
	N_max = N_k.max()
 
  # Create a list of indices of all configurations in kn-indexing.
	mask_kn = zeros([K,N_max], dtype=bool_)
	for k in range(0,K):
	  mask_kn[k,0:N_k[k]] = True

  # Create a list form this mask.
	sample_indices = where(mask_kn)

  # Construct histogram bins
	M = nbins     #number of energy bins
	SMALL = 1.0e-6
	U_min = U_kn[sample_indices].min()
	U_max = U_kn[sample_indices].max()
	U_max += (U_max - U_min) * SMALL  #increment by a bit
	delta_U = (U_max - U_min) / float(M)
	if (bind_width !=None):
		deltat_U = bind_width
		M = int(math.ceil((U_max - U_min) / bin_width))
		print "Using %d bins to achieve bin width of %f" % (M,delta_U)
	else:
		print "Bin width is %f energy units to achieve nbins = %d" % (delta_U, M)
	U_m = U_min + delta_U * (0.5 + arange(0,M,dtype = float64))
	H_m = zeros([M], float64)

  # Assign snapshots to energy bins
	bin_kn = zeros([K,N_max], int32) # bin_kn[k,n] is bin index of U_kn[k,n]
	bin_kn[sample_indices] = array(((U_kn[sample_indices] - U_min) / delta_U), int32)
	H_m = bincount(bin_kn[sample_indices])

  # Compute logs of various quantities
	LOG_ZERO = -1000.0 #replacement for log(0)
	log_H_m = ones([M], float64) * LOG_ZERO
	for m in range(M):
		if (H_m[m] > 0):
			log_H_m[m] = log(H_m[m])
	log_N_k = ones([K], float64) * LOG_ZERO
	for k in range(K):
		if(N_k[k] > 0):
			log_N_k[k] = log(N_k[k])

  # Initialize free energies
	f_k = zeros([K], float64)
	if (initial_f_k != None):
		f_k = initial_f_k.copy()

  # Iterate
	f_k_new = zeros([K], float64)
	max_delta = 0.0
	for iteration in range(maximum_iterations):
  		if (iteration % 5 == 0): print "Iteration %d" % iteration
		
		#Form auxiliary matrices, used in both self-consistent iteration and Newton-Raphson
		# in log space
		log_w_mi = zeros([M,K], float64)
		for m in range(M):
			#denominator = \sum_k N_k exp[f_k - \beta_k U_m]
			exp_arg_k = log_N_k[:] + f_k[:] - beta_k[:]*U_m[m]
			log_w_mi[m,:] = exp_arg_k[:] - logSum(exp_arg_k[:])
		# in real space
		w_mi = zeros([M,K], float64)
		for m in range(M):
			exp_arg_k = f_k[:] - beta_k[:]*U_m[m]
			max_arg = exp_arg_k.max()
			numerators_i = N_k[:] * numpy.exp(exp_arg_k[:] - max_arg)
			w_mi[m,:] = numerators_i[:] / numerators_i.sum()
		
		# Compute new estimates of log weights using self-consistent iteration.
		for i in range(K):
			#Compute new estimate of log weights.
#			f_k_new[i] = f_k[i] + log(N_k[i]) - logSum(log_H_m[:] + log_w_mi[:,i])
			f_k_new[i] = f_k[i] + log_N_k[i] - logSum(log_H_m[:] + log_w_mi[:,i])		
		# Shift to ensure f_k_new[0] = 0.0
		f_k_new -= f_k_new[0]
		# Store difference
		Delta_f_k = f_k_new - f_k
		# update f_k
		f_k = f_k_new.copy()
	
		# If all f_k are zero, terminate
		if all(f_k == 0.0):
			break

		# Terminate when max((f - fold) / f) < relative_tolerance for all nonzero f.
		max_delta = (numpy.abs(Delta_f_k) / (numpy.abs(f_k_new)).max()).max()
		if (iteration % 5 == 0): print "Iteration %8d relative max_delta = %8e" % (iteration, max_delta)
		if isnan(max_delta) or (max_delta < relative_tolerance):
			break

	# if unconverged
	if (iteration == maximum_iterations):
		print "Did not converge in %d iterations (final relative tolerance %e)" % (maximum_iterations, max_delta)
		print "f_k = "
		print f_k
		return f_k
	
	#summary
	print "Converged to relative tolerance of %e (convergence tolerance %e) in %d iterations" % (max_delta, relative_tolerance, iteration)
	print "f_k = "
	print f_k

	# Return the estimate of the dimensionless free energies
	return f_k

def write_free_energies(filename, f_k):
	"""Write free energies to the given file.

	ARGUMENTS
		filename (string) - the name of the file to write to
		f_k (numpy array) - weights to write
	"""

	# get size
	K = f_k.size

	# write out to file
	outfile = open(filename, 'w')
	for k in range(K):
		outfile.write("%f " % f_k[k])
	outfile.write("\n")
	outfile.close()
	
	return

def read_free_energies(filename):
	"""Reads free energies from the given file.
	
	ARGUMENTS
		filename (string) - the name fo the file to write to
	"""
	
	print "Reading free energies from %s..." % filename

	#read file
	infile = open(filename, 'r')
	lines - infile.readlines()
	infile.close()
	
	#split into elements
	elements = list()
	for line in lines:
		elements += line.split()
	print elements

	# parse
	K = len(elements)
	f_k = zeros([K], float64)
	for k in range(K):
		f_k[k] = float(elements[k])

	return f_k

#========================================================================
# MAIN
#========================================================================

#========================================================================
# Read Temperatures
#========================================================================

# Read list of temperatures.
lines0 = read_file(temperature_list_filename)

# Construct list of temperatures (these are strings)
temperatures0 = lines0[0].split()

# Create numpy array of temperatures
num_temps0 = len(temperatures0)   #This used to be K1

# Initialize float vectors to hold temperatures
temperature_vec0 = zeros([num_temps0], float32) # temperature_k[k] is temperature of temperature index k in K

# Convert string list of temperatures to a float vector
	
for n in range(num_temps0):
	temperature_vec0[n] = float(temperatures0[n])

All_Temps = numpy.insert(temperature_vec0, 0, MinT)

# A Loop for each of the three replicates, inserting all temperatures into All_Temps
# Temperatures start at the lowest and move up to the highest, inserting a new one
# at increments specified above
i = 1
j = 0   #temperature_vec0
flag = 0
while (All_Temps[i] <= temperature_vec0[num_temps0-1]):

	if MinTrans < All_Temps[i] < MaxTrans :	increment = increment2
	else : increment = increment1

        while (All_Temps[i-1] + increment < All_Temps[i]):
		if MinTrans < All_Temps[i] < MaxTrans :	increment = increment2
		else : increment = increment1
		All_Temps = numpy.insert(All_Temps, i, All_Temps[i-1] + increment)
		i = i + 1
        i = i + 1
	if i > len(All_Temps)-1: 
		i = i - 1
		break
	

while (All_Temps[i] < MaxT):

	if MinTrans < All_Temps[i] < MaxTrans :	increment = increment2
	else : increment = increment1

	All_Temps = numpy.insert(All_Temps, len(All_Temps), All_Temps[i] + increment)
	i = i + 1

# debuggging
#All_Temps = temperature_vec0
       
total_temps = All_Temps.size

# Compute inverse temperatures
beta_vec = (kB * All_Temps)**(-1)


#=======================================================================
# Read potential energies
#=======================================================================
print "Reading potential energies..."
lines = read_file(potential_energies_filename)
num_energies = len(lines)
U_kt1 = zeros([num_temps0,num_energies], float32) # U_kt[k,t] is the potential energy (kJ/mol) for snapshot t of temperature index k
print "%d lines read, processing %d snapshots" % (len(lines), num_energies)
for t in range(num_energies):
	# Get line containing the energies for snapshot t
	line = lines[t]
	# Extract energy values from text
	elements = line.split()
	for k in range(num_temps0):
		U_kt1[k,t] = float(elements[k])

# Initialize expanded U_kt for new temperatures with zeros
U_kt = zeros([total_temps,num_energies], float32)
N_k = zeros([total_temps], int32)

# Copy potential energies into columns where old and new temperatures match
y = 0
for j in range(total_temps):
	if All_Temps[j] == temperature_vec0[y]:
		U_kt[j,0:num_energies] = U_kt1[y,0:num_energies]
		N_k[j] = num_energies
		y = y +1
		if y > (num_temps0-1) : y = 0

#U_kn = U_kt.copy() # In case the data needed to be "uncorrelated", there used to be a more sophisticated transfer from U_kt to U_kn, the latter is used from here on
#N_max = T



#========================================================================
# Compute reduced potential energies
#========================================================================

print "Computing reduced potential energies..."
u_kln = zeros([total_temps,total_temps,num_energies], float32) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
for k in range(total_temps):
      for l in range(total_temps):
            u_kln[k,l,0:N_k[k]] = beta_vec[l] * U_kt[k,0:N_k[k]]

#========================================================================
# Use WHAM function call to compute initial guess of f_k
# Initialize MBAR
#========================================================================

free_energies_filename = os.path.join(out_directory, 'f_k.out')
# Compute initial guess using histogram WHAM
#print "Using WHAM to generate histogram-based initial guess of dimensionless free energies f_k..."
#f_k = histogram_wham(beta_vec, U_kt, N_k)

# Initialize MBAR with Newton-Raphson
print "Initializing MBAR......"
mbar = pymbar.MBAR(u_kln, N_k, method = 'adaptive', verbose = True, initialize='BAR')
#mbar = pymbar.MBAR(u_kln, N_k, method = 'self-consistent-iteration', initialize='BAR', verbose = True)
# Write free energies
write_free_energies(free_energies_filename, mbar.f_k)

print "MBAR routines completed, Preparing Data to Compute Expectations..."
print "."
print "."
print "."

#========================================================================
# Prep E_kt
#========================================================================

print "Reading total energies for each temp..."
lines = read_file(total_energies_filename)
num_totE = len(lines)
Ei_kt = zeros([num_temps0,num_totE], float32)
print "Found %d snapshots, processing..." % (num_totE)
for t in range(num_totE):
	line = lines[t]
	elements = line.split()
	for k in range(num_temps0):
		Ei_kt[k,t] = float(elements[k])

E_kt = zeros([total_temps,num_totE], float32)

y = 0
for j in range(total_temps):
	if All_Temps[j] == temperature_vec0[y]:
		E_kt[j,0:num_totE] = Ei_kt[y,0:num_totE]
		y = y+1
		if y > (num_temps0-1) : y = 0

#========================================================================
# Prep E2_kt
#========================================================================

E2_kt = zeros([total_temps,num_totE], float32)
for k in range(total_temps):
      for l in range(total_temps):
            E2_kt[k,0:N_k[k]] = E_kt[k,0:N_k[k]] * E_kt[k,0:N_k[k]];

#=========================================================================
# Compute Expectations for E_kt and E2_kt as E_expect and E2_expect
#=========================================================================

# Define path and file names
E_filename = os.path.join(out_directory, 'E_expect.out')
dE_filename = os.path.join(out_directory, 'dE_expect.out')
E2_filename = os.path.join(out_directory, 'E2_expect.out')
dE2_filename = os.path.join(out_directory,'dE2_expect.out')
K_filename = os.path.join(out_directory, 'temperatures.out')

print "Computing Expectations for E and E^2..."
(E_expect, dE_expect) = mbar.computeExpectations(E_kt)
(E2_expect,dE2_expect) = mbar.computeExpectations(E2_kt)

# Output expectations to file
print "Outputting Files for E and E^2..."
write_free_energies(E_filename, E_expect)
write_free_energies(dE_filename, dE_expect)
write_free_energies(E2_filename, E2_expect)
write_free_energies(dE2_filename, dE2_expect)
write_free_energies(K_filename, All_Temps)

#==========================================================================
# Read in contact info
#==========================================================================

print "Reading contact ratios for each temp..."
lines = read_file(contacts_filename)
num_contacts = len(lines)
C1_kt = zeros([num_temps0,num_contacts], float32)
print "Found %d snapshots, processing..." % (num_contacts)
for t in range(num_contacts):
	line = lines[t]
	elements = line.split()
	for k in range(num_temps0):
		C1_kt[k,t] = float(elements[k])

C_kt = zeros([total_temps,num_contacts], float32)
y = 0
for j in range(total_temps):
	if All_Temps[j] == temperature_vec0[y]:
		C_kt[j,0:num_contacts] = C1_kt[y,0:num_contacts]
		y = y+1
		if y > (num_temps0-1) : y = 0

#===========================================================================
# Compute expectations for C_kt, or contact ratios
#===========================================================================
C_filename = os.path.join(out_directory, 'contacts.out')
d_C_filename = os.path.join(out_directory, 'd_contacts.out')

print "Computing Expectations for Contact Ratios..."
(C_expect, d_C_expect) = mbar.computeExpectations(C_kt)

print "Writing contact expectations to file..."
write_free_energies(C_filename, C_expect)
write_free_energies(d_C_filename, d_C_expect)

print "Done!"
