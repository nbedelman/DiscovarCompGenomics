""""
performs simulation of specified tree for each chromosome, making a directory to store it in based on the input

use in the following way

python runsim_scenario4.py int modelname
"""

import subprocess
import msprime
import math
import numpy as np
import sys


##subprocess.call(['mkdir', sys.argv[1]])
filename = "%s/%s" %(sys.argv[1], sys.argv[2])
mutation_rate = 3e-9
recombination_rate = 0#3e-9
#length = 5e8
sample_size = 1
max_chrom = 1000
pop_size = 1000000
mig_rate = 0.2
theta=2*2000000*3e-9*1e3
length =1e3
lengthStore='1e3'
#print theta

def run_sim(chrom, length, mutation_rate, recombination_rate, sample_size, pop_size, filename, mig_rate):
	# M is the overall symmetric migration rate, and d is the number
	# of demes.
	M = 0.0
	d = 4
    	# We rescale m into per-generation values for msprime.
    	m = M / (4 * (d - 1))
   	
        # Population Size
        N_1 = pop_size
        N_2 = pop_size
        N_3 = pop_size
	N_4 = pop_size

	# Allocate the initial sample. Because we are interested in the
    	population_configurations = [
        	msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_1),
       		msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_2),
        	msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_3),
		msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_4)]


   	
	# Now we set up the migration matrix.
  	migration_matrix = [
        	[0, 0, 0, 0],
        	[0, 0, 0, 0],
        	[0, 0, 0, 0],
		[0, 0, 0, 0]]
    
    #Times are given in generations
    	T_3_2 = 6000000
    	T_2_1 = 3000000
	T_admix = 2000000
	T_mrca = 20000000
    
    	demographic_events = [


	# One-way flow from 1 into 3
	# 1 splits from 2, i.e. 1 and 2 merge
		msprime.MassMigration(
			time = T_admix, source = 0, destination = 2, proportion = 0.3),
		msprime.MassMigration(
			time = T_2_1, source = 0, destination = 1, proportion = 1.0),
	       # msprime.MigrationRateChange(
        	#    	time=T_2_1, rate=0, matrix_index=(0, 1)),
                #msprime.MigrationRateChange(
                  #      time=T_2_1, rate=0, matrix_index=(1, 0)),
		#msprime.PopulationParametersChange(
		#	time = T_2_1, initial_size=1, growth_rate = 0, population_id = 0),

	# There's one-way gene flow from population 1 into 3.

	# 2 splits from 3, i.e. 2 and 3 merge
		msprime.MassMigration(
			time = T_3_2, source = 1, destination = 2, proportion = 1.0),
	#3 and 4 merge:
		msprime.MassMigration(
			time = T_mrca, source = 2, destination = 3, proportion = 1.0)
	
               # msprime.MigrationRateChange(
               #         time=T_3_2, rate=0, matrix_index=(1, 2)),
                #msprime.MigrationRateChange(
                #        time=T_3_2, rate=0, matrix_index=(2, 1)),
		#msprime.PopulationParametersChange(
		#	time = T_3_2, initial_size=1, growth_rate = 0, population_id = 1),			
			
	]			
			

#	# Use the demography debugger to print out the demographic history
#    	# that we have just described.
   	dp = msprime.DemographyDebugger(
    		Ne=N_1,
	       	population_configurations=population_configurations,
        	migration_matrix=migration_matrix,
        	demographic_events=demographic_events)
    	dp.print_history()

	dp = msprime.simulate(
		Ne=N_1,
		population_configurations = population_configurations,
		migration_matrix = migration_matrix,
		demographic_events = demographic_events,
		mutation_rate = mutation_rate, 
		length = length,
		recombination_rate = recombination_rate)
#	shape = dp.get_num_mutations(), dp.get_sample_size()
#	A = np.empty(shape, dtype = "u1")	

	#print dp.first().newick()
	output=''
	runningTot=0
	for index,tr in enumerate(dp.trees()):
		lengthQ=-int(np.round(tr.get_interval()[0]))+int(np.round(tr.get_interval()[1]))
		if lengthQ>0:
			runningTot+=lengthQ
			output=output+'['+str(lengthQ)+']'+str(tr.newick())+'\n'
	#print runningTot

	#generate the output eigenstrat files
	seqgeninput=open('%s_%s_%s.txt' %(filename,chrom,lengthStore),'w')
	seqgeninput.write(output)
	seqgeninput.close()
#	geno = open('%s_chr%s.geno' %(filename, chrom), 'w')
#	snp = open('%s_chr%s.snp' %(filename, chrom), 'w')
#	for variant in dp.variants(as_bytes=True):	
#		geno.write('%s\n' %variant.genotypes)
#		snp.write('%s	%s	%s	%s	0	1\n' %(variant.index, chrom, variant.position/float(length), variant.position))
#	geno.close()
#	np.close()

#	totalsamplesize = dp.get_sample_size()
#	ind = open('%s_chr%s.ind' %(filename, chrom), 'w')
#	pop = -1
#	for row in range(0, totalsamplesize):
#		if row % sample_size == 0:
#			pop += 1			
#		ind.write('%s	U	%s\n' %(row, pop))
#	ind.close()

# run simulation for each chromosome
for chrom in range(1,max_chrom+1):
	run_sim(chrom, length, mutation_rate, recombination_rate, sample_size, pop_size, filename, mig_rate)
