import subprocess
import msprime
import math
import numpy as np
import sys

#subprocess.call(['mkdir', sys.argv[1]])
#subprocess.call(['mkdir', str(sys.argv[1]+'/'+str(sys.argv[2])])
#filename = "%s/%s" %(sys.argv[1], sys.argv[2])
#mutation_rate = 3e-9
#recombination_rate = 0
#sample_size = 1
#max_chrom = 5000
#admix_prop = float(sys.argv[2])#0.3
##theta=2*2000000*3e-9*1e3
#length =1e3
#lengthStore='1e3'
#pop_size = int(sys.argv[1])#1000000
##print theta

nlist=[250000]#[500000,1000000,2000000]
proplist=[0.2,0.4,0.6,0.8]

def run_sim(chrom, length, mutation_rate, recombination_rate, sample_size, pop_size, filename, admix_prop):
	d = 4	
        # Population Size
        N_1 = pop_size
        N_2 = pop_size
        N_3 = pop_size
	N_4 = pop_size

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
	T_mrca = 10000000
    
    	demographic_events = [


	# One-way flow from 2 into 3
	# 1 splits from 2, i.e. 1 and 2 merge
		msprime.MassMigration(
			time = T_admix, source = 1, destination = 2, proportion = admix_prop),
		msprime.MassMigration(
			time = T_2_1, source = 0, destination = 1, proportion = 1.0),
	       # msprime.MigrationRateChange(
        	#    	time=T_2_1, rate=0, matrix_index=(0, 1)),
                #msprime.MigrationRateChange(
                  #      time=T_2_1, rate=0, matrix_index=(1, 0)),
		msprime.PopulationParametersChange(
			time = T_2_1, initial_size=1, growth_rate = 0, population_id = 0),

	# 2 splits from 3, i.e. 2 and 3 merge
		msprime.MassMigration(
			time = T_3_2, source = 1, destination = 2, proportion = 1.0),
	#3 and 4 merge:
		msprime.MassMigration(
			time = T_mrca, source = 2, destination = 3, proportion = 1.0)
	]			
			

	dp = msprime.simulate(
		Ne=N_1,
		population_configurations = population_configurations,
		migration_matrix = migration_matrix,
		demographic_events = demographic_events,
		mutation_rate = mutation_rate, 
		length = length,
		recombination_rate = recombination_rate)
	output=''
	runningTot=0
	for index,tr in enumerate(dp.trees()):
		lengthQ=-int(np.round(tr.get_interval()[0]))+int(np.round(tr.get_interval()[1]))
		if lengthQ>0:
			runningTot+=lengthQ
			output=output+'['+str(lengthQ)+']'+str(tr.newick())+'\n'
	seqgeninput=open('%s_%s.txt' %(filename,chrom),'w')
	seqgeninput.write(output)
	seqgeninput.close()

for n in nlist:
	for p in proplist:
		subprocess.call(['mkdir', str(n)])
		subprocess.call(['mkdir', str(n)+'/'+str(p)])
		filename = "%s/%s/%s_%s" %(str(n), str(p),str(n),str(p))
		mutation_rate = 3e-9
		recombination_rate = 0
		sample_size = 1
		max_chrom = 5000
		admix_prop = float(p)#0.3
		#theta=2*2000000*3e-9*1e3
		
		length =1e3
		lengthStore='1e3'
		pop_size = int(n)#1000000
		for chrom in range(1,max_chrom+1):
			run_sim(chrom, length, mutation_rate, recombination_rate, sample_size, pop_size, filename, admix_prop)



# run simulation for each chromosome
#for chrom in range(1,max_chrom+1):
#	run_sim(chrom, length, mutation_rate, recombination_rate, sample_size, pop_size, filename, mig_rate)
