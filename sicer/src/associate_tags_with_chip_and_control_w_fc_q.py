import re, os, sys, shutil
from math import *
from string import *
from decimal import Decimal

import numpy as np
import multiprocessing as mp
from functools import partial

from sicer.lib import GenomeData;
from sicer.lib import associate_tags_with_regions

import scipy
import scipy.stats


def associate_tag_count_to_regions(args,scaling_factor, control_library_size, genomesize, chrom):
	island_file = args.treatment_file.replace('.bed','')+'_'+chrom+'_graph.npy'
	treatment_file = args.treatment_file.replace('.bed','')+'_'+chrom+'.npy'
	control_file = args.control_file.replace('.bed','')+'_'+chrom+'.npy'

	island_list = np.load(island_file)
	island_start_list = []
	island_end_list = []

	for island in island_list:
		island_start_list.append(island[1])
		island_end_list.append(island[2])

	total_chip_count = 0
	island_chip_readcount_list = [0]*len(island_list)
	treatment_reads = np.load(treatment_file)
	for read in treatment_reads:
		position = associate_tags_with_regions.tag_position(read,args.fragment_size)
		index = associate_tags_with_regions.find_readcount_on_islands(island_start_list,island_end_list,position)
		if index >= 0: #if the read is found in a valid island
			island_chip_readcount_list[index] += 1
			total_chip_count+=1

	total_control_count = 0
	island_control_readcount_list = [0]*len(island_list)
	control_reads = np.load(control_file)
	for read in control_reads:
		position = associate_tags_with_regions.tag_position(read, args.fragment_size)
		index = associate_tags_with_regions.find_readcount_on_islands(island_start_list,island_end_list,position)
		if index >=0:
			island_control_readcount_list[index]+=1
			total_control_count+=1

	output_lines = []
	pvalue_list = []
	for index in range(0,len(island_list)):
		island = island_list[index]
		observation_count = island_chip_readcount_list[index]
		control_count = island_control_readcount_list[index]
		average=0
		if(control_count>0):
			average = control_count * scaling_factor
		else:
			length = island[2]-island[1]+1
			average = length * control_library_size*1.0/genomesize
			average = min(0.25, average)* scaling_factor;
		fc = float(observation_count)/float(average)
		if (observation_count>average):
			pvalue = scipy.stats.poisson.sf(observation_count, average)
		else:
			pvalue = 1

		pvalue_list.append(pvalue)
		output_line = [island[0],island[1],island[2],observation_count,control_count, pvalue,fc]
		output_lines.append(output_line)

	np_output_lines = np.array(output_lines,dtype=object)
	file_name = args.treatment_file.replace('.bed','')+'_'+chrom+'_'+'island_summary.npy'
	np.save(file_name,np_output_lines)
	pvalue_save_name = chrom+'_pvalue.npy'
	np.save(pvalue_save_name,pvalue_list)
	return pvalue_save_name


def main(args,chip_library_size, control_library_size):

	chroms = GenomeData.species_chroms[args.species];
	genomesize = sum (GenomeData.species_chrom_lengths[args.species].values());
	genomesize = args.effective_genome_fraction* genomesize;

	print ("chip library size:  ", chip_library_size)
	print ("control library size:  ", control_library_size)

	totalchip = 0;
	totalcontrol = 0;
	scaling_factor = chip_library_size*1.0/control_library_size

	island_chip_readcount = {};
	island_control_readcount = {};

	#Use multiprocessing to associate each read with an island
	pool = mp.Pool(processes = min(mp.cpu_count(),len(chroms)))
	associate_tag_count_to_regions_partial = partial(associate_tag_count_to_regions,args,scaling_factor,control_library_size,genomesize)
	p_value_files = pool.map(associate_tag_count_to_regions_partial, chroms)
	pool.close()

	#Get the list of p-value from each parallel processes and concatenate them into one list of all p-values
	p_value_list = np.array([])
	for p_value_file in p_value_files:
		chrom_p_value_list=np.load(p_value_file)
		p_value_list = np.concatenate([p_value_list,chrom_p_value_list])
		os.remove(p_value_file)
	p_value_rank_array = scipy.stats.rankdata(p_value_list)
	total_num_of_pvalue = len(p_value_list)


	index = 0
	file_name = args.treatment_file.replace('.bed','')
	outfile_path = os.path.join(args.output_directory,(file_name+'-W'+str(args.window_size)+'-G'+str(args.gap_size)+'-islands-summary'))
	with open(outfile_path,'w') as outfile:
		for chrom in chroms:
			island_file_name = file_name+'_'+chrom+'_'+'island_summary.npy'
			island = np.load(island_file_name)
			modified_island = []
			for i,line in enumerate(island):
				totalchip+=int(line[3])
				totalcontrol+=int(line[4])
				alpha_stat = p_value_list[index] * total_num_of_pvalue/p_value_rank_array[index];
				if alpha_stat > 1:
					alpha_stat = 1;

				line=line.tolist()
				line.append(alpha_stat)
				outputline=(line[0]+'\t'+str(line[1])+'\t'+str(line[2])+'\t'+str(line[3])+'\t'+str(line[4])+'\t'+
							str(line[5])+'\t'+str(line[6])+'\t'+str(line[7])+'\n')
				outfile.write(outputline)
				modified_island.append(tuple(line))
				index+=1

			np_modified_island=np.array(modified_island,dtype=object)
			np.save(island_file_name,np_modified_island)

	print("Total number of chip reads on islands is: ", totalchip)
	print("Total number of control reads on islands is: ",totalcontrol)
if __name__ == "__main__":
	main(sys.argv)
