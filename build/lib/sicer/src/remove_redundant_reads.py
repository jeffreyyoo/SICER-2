
import re, os, sys, shutil
from math import *
from string import *

import subprocess
import multiprocessing as mp
from functools import partial
import numpy as np

from sicer.lib import GenomeData

'''Filters redundant reads according to the cutoff value by taking a sorted list and comparing adjacent reads'''
def remove_redundant_1chrom_single_strand_sorted(reads,cutoff):
	current_start = 0
	current_end = 0
	current_count = 1
	total = 0
	retained = 0
	filtered_reads=[]

	for i,read in enumerate(reads):
		total += 1
		start = read[1]
		end = read[2]
		if start != current_start:
			retained += 1
			current_start = start
			current_end = end
			current_count = 1
			filtered_reads.append(read)
		elif end != current_end:
			retained += 1
			current_start = start
			current_end = end
			current_count = 1
			filtered_reads.append(read)
		else:
			current_count += 1
			if current_count <= cutoff:
				retained += 1
				filtered_reads.append(read)

	return (total, retained, filtered_reads)

'''Separates reads by positive and negative strands before filtering redudant reads.
	Saves the filtered reads as a numpy binary file in temporary directory created in run_SICER.
	This is because python cannot pass extremely large objects between parallel processes'''
def strand_broken_remove(chrom, cutoff, file, chrom_reads):
	#Use of multiprocessing means print statements will be out of order. Use print_return to hold them until the end
	print_return = ""
	plus_reads=[]
	minus_reads=[]
	for read in chrom_reads:
		if(read[5]=='+'):
			plus_reads.append(read)
		elif(read[5]=='-'):
			minus_reads.append(read)
	if(not plus_reads):
		sys.stderr.write(chrom+" + reads do not exist in "+file+"\n")
	if(not minus_reads):
		sys.stderr.write(chrom+" - reads do not exist in "+file+"\n")

	plus_reads = sorted(plus_reads,key = lambda x: (x[1],x[2]))
	minus_reads = sorted(minus_reads,key = lambda x: (x[1],x[2]))
	(p_total, p_retained,filtered_plus_reads) = remove_redundant_1chrom_single_strand_sorted(plus_reads, cutoff)
	(m_total, m_retained,filtered_minus_reads) = remove_redundant_1chrom_single_strand_sorted(minus_reads, cutoff)

	print_return += (chrom+ "\tPlus reads: "+str(p_total)+"\tRetained plus reads: "+str(p_retained)+";\tMinus reads: "
					+str(m_total)+"\tRetained minus reads: "+str(m_retained))

	filtered_output = filtered_plus_reads + filtered_minus_reads
	np_filtered_output = np.array(filtered_output,dtype=object)
	name_for_save = file+"_"+chrom+".npy"
	np.save(name_for_save,np_filtered_output)
	total_retained = p_retained + m_retained

	return (print_return,total_retained)


'''Opens the given file and greps all the reads of the given chromsome and stores them in a list of tuples
	Tuple format: (chrom,start,end,name,score,strand)'''
def match_by_chrom(file, chrom):
	match = chrom + "[[:space:]]"
	matched_reads = subprocess.Popen('grep %s %s' % (match,file), stdout=subprocess.PIPE,shell=True)
	chrom_reads= str(matched_reads.stdout.read(),'utf-8').splitlines()  #generates a list of each reads, which are represented by a string value

	for i,reads in enumerate(chrom_reads):
		reads=re.split('\t', reads)
		reads[1]=int(reads[1])
		reads[2]=int(reads[2])
		chrom_reads[i]=tuple(reads)

	return chrom_reads


'''Function designed for handling multiprocessing. Separates all reads by chromosome
	and then filters redudant reads'''
def find_and_filter_reads (path_to_file, cutoff, chrom):
	file_name = os.path.basename(path_to_file)
	file_name = file_name.replace('.bed','')

	chrom_reads = match_by_chrom(path_to_file, chrom) #Separates all reads by chromosome
	return strand_broken_remove(chrom, cutoff, file_name, chrom_reads)




'''path_to_file: complete path to the .bed file that needs to processed for redudant reads'''
def main(args,path_to_file):
	chroms = GenomeData.species_chroms[args.species];  #list of chromsomes of the given species
	cutoff = args.redundancy_threshold

	#Use multiprocessing module to run parallel processes for each chromosome
	pool = mp.Pool()
	find_and_filter_reads_partial= partial(find_and_filter_reads,path_to_file, cutoff)
	filtered_result = pool.map(find_and_filter_reads_partial, chroms)
	pool.close()


	total_read_count = 0
	for result in filtered_result:
		print(result[0])
		total_read_count+=result[1]

	return total_read_count


if __name__ == "__main__":
	main(sys.argv,"")
