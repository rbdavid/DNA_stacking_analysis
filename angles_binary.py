#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import MDAnalysis
import sys
import os
import matplotlib.pyplot as plt

traj_file ='%s' %(sys.argv[1])

# ----------------------------------------
# VARIABLE DECLARATION

base1 = 1
nbases = 15
#nbases = 3 

#Nsteps = 150000		# check length of the energy file; if not 150000 lines, then need to alter Nsteps value so that angle values will match up
#Nsteps = 149996
#equilib_step = 37500		# we have chosen 75 ns to be the equilib time; 75ns = 37500 frames; if energy values do not match with angle values, then equilib_step needs to be altered as well...
#equilib_step = 37496

#production = Nsteps - equilib_step

# SUBROUTINES/DEFINITIONS:

arccosine = np.arccos
dotproduct = np.dot
pi = np.pi

ldtxt = np.loadtxt
zeros = np.zeros

# ----------------------------------------
# DICTIONARY DECLARATION 

normals = {}			# create the normals dictionary for future use
total_binaries = {}		# create the total_binaries dictionary for future use
get_norm = normals.get
get_tb = total_binaries.get

# ----------------------------------------
# PLOTTING SUBROUTINES

def plotting(xdata, ydata, base):
	plt.plot(xdata, ydata, 'rx')
	plt.title('Stacking behavior of base %s over the trajectory' %(base))
	plt.xlabel('Simulation time (ns)')
	plt.ylabel('Stacking metric')
	plt.xlim((0,300))
	plt.grid( b=True, which='major', axis='both', color='k', linestyle='-')
	plt.savefig('stacking_binary.%s.png' %(base))
	plt.close()

def vdw_hist(data, base_a, base_b):
	events, edges, patches = plt.hist(data, bins = 100, histtype = 'bar')
	plt.title('Distribution of vdW Energies - Base Pair %s-%s' %(base_a, base_b))
	plt.xlabel('vdW Energy ($kcal\ mol^{-1}$)')
	plt.xlim((-8,0))
	plt.ylabel('Frequency')
	plt.savefig('energy.%s.%s.png' %(base_a, base_b))
	nf = open('energy.%s.%s.dat' %(base_a, base_b), 'w')
	for i in range(len(events)):
		nf.write(' %10.1f    %10.4f\n' %(events[i], edges[i]))
	nf.close()
	plt.close()
	events = []
	edges = []
	patches = []

def angle_hist(data, base_a, base_b):
	events, edges, patches = plt.hist(data, bins = 100, histtype = 'bar')
	plt.title('Distribution of Angles btw Base Pair %s-%s' %(base_a, base_b))
	plt.xlabel('Angle (Degrees)')
	plt.ylabel('Frequency')
	plt.savefig('angle.%s.%s.png' %(base_a, base_b))
	nf = open('angle.%s.%s.dat' %(base_a, base_b), 'w')
	for i in range(len(events)):
		nf.write(' %10.1f    %10.4f\n' %(events[i], edges[i]))
	nf.close()
	plt.close()
	events = []
	edges = []
	patches = []

def energy_angle_hist(xdata, ydata, base_a, base_b):
	counts, xedges, yedges, image = plt.hist2d(xdata, ydata, bins = 100)
	cb1 = plt.colorbar()
	cb1.set_label('Frequency')
	plt.title('Distribution of Base Pair interactions - %s-%s' %(base_a, base_b))
	plt.xlabel('Angle (Degrees)')
	plt.ylabel('vdW Energy ($kcal\ mol^{-1}$)')
	plt.ylim((-6,0.5))
	plt.savefig('vdw_angle.%s.%s.png' %(base_a, base_b))
	plt.close()
	counts = []
	xedges = []
	yedges = []
	image = []

# MAIN PROGRAM:

# ----------------------------------------
# ATOM SELECTION - load the trajectory and select the desired nucleotide atoms to be analyzed later on

u = MDAnalysis.Universe('../nucleic_ions.pdb', traj_file, delta=2.0)    # load in trajectory file

Nsteps = len(u.trajectory)
equilib_step = 37500				# first 75 ns are not to be included in total stacking metric
production = Nsteps - equilib_step

nucleic = u.selectAtoms('resid 1:15')           # atom selections for nucleic chain

a1 = nucleic.selectAtoms('resid 1')             # residue 1 has different atom IDs for the base atoms
a1_base = a1.atoms[10:24]                       # atom selections

bases = []                                      # make a list of the 15 bases filled with atoms
bases.append(a1_base)                           # add base 1 into list

for residue in nucleic.residues[1:15]:          # collect the other bases into list
	residue_base = []
	residue_base = residue.atoms[12:26]
	bases.append(residue_base)

# ----------------------------------------
# DICTIONARY DEVELOPMENT - Develop the normals and total binary dictionary which contain the data for each base

while base1 <= nbases:
	normals['normal.%s' %(base1)] = get_norm('normal.%s' %(base1), np.zeros((Nsteps, 3)))
	total_binaries['base.%s' %(base1)] = get_tb('base.%s' %(base1), np.zeros(Nsteps))
	base1 += 1

# ----------------------------------------
# SIMULATION TIME - calculate the array that contains the simulation time in ns units

time = np.zeros(Nsteps)

for i in range(Nsteps):
	time[i] = i*0.002		# time units: ns

# ----------------------------------------
# NORMAL ANALYSIS for each base - loops through all bases and all timesteps of the trajectory; calculates the normal vector of the base atoms

base1 = 1

while (base1 <= nbases):
	
	for ts in u.trajectory:
		Princ_axes = []
		Princ_axes = bases[base1 - 1].principalAxes()
		normals['normal.%s' %(base1)][ts.frame - 1] = Princ_axes[2]		# ts.frame index starts at 1; add normal to dictionary with index starting at 0
	
	base1 += 1

# ----------------------------------------
# BASE PAIR ANALYSIS - loops through all base pairs (w/out duplicates) and performs the angle analysis as well as the binary analysis

base1 = 1		# reset the base index to start at 1

while (base1 <= nbases):		# while loops to perform the base-pair analysis while avoiding performing the same analysis twice
	base2 = base1 + 1
	while (base2 <= nbases):
		
		os.mkdir('base%s_base%s' %(base1, base2))		# makes and moves into a directory for the base pair
		os.chdir('base%s_base%s' %(base1, base2))

		energyfile = '../../nonbond_energy/base%s_base%s/base%s_base%s.energies.dat' %(base1, base2, base1, base2)
		
		energies = ldtxt(energyfile)		# load in the energy file to a numpy array
		
		vdw_energies = energies[:,2]

		binary = zeros(Nsteps)
		nf = open('binary.%s.%s.dat' %(base1, base2), 'w')		# write the base pair data to a file; make sure to be writing this in a base pair directory

		# angle and binary analysis for base pair; 
		for i in range(Nsteps):
			
			angle = 0.
			angle = arccosine(dotproduct(normals['normal.%s' %(base1)][i], normals['normal.%s' %(base2)][i]))
			angle = angle*(180./pi)
			if angle > 90.:
				angle = 180. - angle
			
			if vdw_energies[i] <= -3.5 and angle <= 30.:	# cutoff: -3.5 kcal mol^-1 and 30 degrees
				binary[i] = 1.					# assumed else binary[i] = 0.
			
			nf.write(' %10.3f    %10.5f    %10.5f    %10.1f\n' %(time[i], vdw_energies[i], angle, binary[i])) # check time values
			
			total_binaries['base.%s' %(base1)][i] = total_binaries['base.%s' %(base1)][i] + binary[i]
			total_binaries['base.%s' %(base2)][i] = total_binaries['base.%s' %(base2)][i] + binary[i]

		nf.close()
		
		angles = []
		energies = []
		vdw_energies = []

		os.chdir('..')

		base2 += 1

	base1 += 1

# ----------------------------------------
# TOTAL BINARY METRIC ANALYSIS - writing to file and plotting 
# print out (also plot) the total binary data to an indivual file for each individual base

base1 = 1		# reset the base index to start at 1

os.mkdir('total_binaries')
os.chdir('total_binaries')

while (base1 <= nbases):
	
	os.mkdir('base%s' %(base1))
	os.chdir('base%s' %(base1))

	nf = open('binary.%s.dat' %(base1), 'w')

	for i in range(Nsteps):
		
		nf.write(' %10.3f    %10.1f\n' %(time[i], total_binaries['base.%s' %(base1)][i]))		# check time values

	nf.close()
	
	counts = 0
	
	for i in range(equilib_step, Nsteps):
		if total_binaries['base.%s' %(base1)][i] > 0.:
			counts +=1
	
	prob = 0.
	prob = (float(counts)/production)*100.

	nf = open('stacking.%s.dat' %(base1), 'w')
	nf.write('counts: %10.1f out of %10.1f time steps \n Probability of stacking = %10.4f ' %(counts, production, prob))
	nf.close()

	plotting(time[:], total_binaries['base.%s' %(base1)][:], base1)
	
	os.chdir('..')

	base1 += 1

# ----------------------------------------
# BASE PAIR PLOTTING - making histogram plots for vdW energy distributions, angle distributions, and 2d hist of vdw vs angle distributions
# Also printint out a file that contains the count of timesteps where the base pair are stacked

os.chdir('..')

base1 = 1

while (base1 <= nbases):                # while loops to perform the base-pair analysis while avoiding performing the same analysis twice
        base2 = base1 + 1
        while (base2 <= nbases):
		os.chdir('base%s_base%s' %(base1, base2))
		
		infile = 'binary.%s.%s.dat' %(base1, base2)
		data = ldtxt(infile)		# data[0] = time, data[1] = vdW energies, data[2] = angle, data[3] = base pair binary metric
		
		vdw_hist(data[equilib_step:,1], base1, base2)

		angle_hist(data[equilib_step:,2], base1, base2)

		energy_angle_hist(data[equilib_step:,2], data[equilib_step:,1], base1, base2)

		nf = open('stacking.%s.%s.dat' %(base1, base2), 'w')
		bp_counts = sum(data[equilib_step:,3])
		nf.write('counts for base pair %s-%s: %10.1f' %(base1, base2, bp_counts))
		nf.close()

		data = []

		os.chdir('..')

		base2 += 1
	
	base1 += 1

# ----------------------------------------
# END
