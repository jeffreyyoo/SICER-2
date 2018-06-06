#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
#
# Disclaimer
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010


import re, os, sys, shutil
from math import *
from string import *

import multiprocessing as mp
from functools import partial
import numpy as np

from lib import GenomeData


def write(item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside.
	item is a BED3 object
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\n";
	out.write(outline);

# #Opens the given file and greps all the reads of the given chromsome and stores them in a list of tuples
# #Tuple format: (chrom,start,end,score/count)
# def match_by_chrom(file, chrom):
# 	match = chrom + "[[:space:]]"
# 	matched_reads = subprocess.Popen('grep %s %s' % (match,file), stdout=subprocess.PIPE,shell=True)
# 	chrom_reads= str(matched_reads.stdout.read(),'utf-8').splitlines()
#
# 	for i,reads in enumerate(chrom_reads):
# 		reads=re.split('\t', reads)
# 		reads[1]=int(reads[1])
# 		reads[2]=int(reads[2])
# 		reads[3]=float(reads[3])
# 		chrom_reads[i]=tuple(reads)
#
# 	return chrom_reads

#Function designed for handling multiprocessing. Executes the redundancy removal algorithm
#for each independent chromosome
def find_union_islands (no_control, temp_dir_1, temp_dir_2, chrom ):
	file_name_1 = temp_dir_1 + '_'+chrom
	file_name_2 = temp_dir_2 + '_'+chrom
	if(no_control == True):
		file_name_1 += '_graph.npy'
		file_name_2 += '_graph.npy'
	else:
		file_name_1 += '_island_summary.npy'
		file_name_2 += '_island_summary.npy'

	island_list_1 = np.load(file_name_1)
	island_list_2 = np.load(file_name_2)
	if(len(island_list_1)==0):
		island_list = island_list_2
	elif(len(island_list_2)==0):
		island_list = island_list_1
	else:
		island_list = np.concatenate((island_list_1, island_list_2))

	union_island_list = []
	if(len(island_list)>0):
		island_list = island_list[np.argsort(island_list[:,1])]
		current = island_list[0]
		i = 1
		while i < len(island_list):
			compare = island_list[i]
			assert current[1] <= compare[2]
			if compare[1] > current[2]:
				union_island_list.append(current)
				current = compare
				i += 1
			else:
				current[2] = max(current[2], compare[2])
				i += 1
		union_island_list.append(current)
	np_union_island_list = np.array(union_island_list, dtype=object)
	np.save(chrom+'_union_output.npy',np_union_island_list)

def main(args,temp_dir_1,temp_dir_2):
	chroms = GenomeData.species_chroms[args.species]
	# path_to_file_1 = args.output_directory
	# path_to_file_2 = args.output_directory
	# file_name_1 = args.treatment_file[0].replace('.bed','')
	# file_name_2 = args.treatment_file[1].replace('.bed','')
	# if(len(args.control_file)==0):	#If there is no control libraries
	# 	path_to_file_1 += ('/'+file_name_1+'-W'+str(args.window_size)
	# 					+'-G'+str(args.gap_size)+'.scoreisland')
	# 	path_to_file_2 += '/'+file_name_2+'-W'+str(args.window_size)
	# 					+'-G'+str(args.gap_size)+'.scoreisland')
	# else:
	# 	path_to_file_1 += ('/'+file_name_1+'-W'+str(args.window_size)+'-G'
	# 					+str(args.gap_size)+'-FDR'+str(args.false_discovery_rate)+'-island.bed')
	# 	path_to_file_2 += ('/'+file_name_2+'-W'+str(args.window_size)+'-G'
	# 					+str(args.gap_size)+'-FDR'+str(args.false_discovery_rate)+'-island.bed')

	#Partially fill out the full directory of the files we want to access
	temp_dir_1 += '/'+args.treatment_file[0].replace('.bed','')
	temp_dir_2 += '/'+args.treatment_file[1].replace('.bed','')

	no_control=False #Default value
	if(len(args.control_file)==0): #if there are no control libraries
		no_control = True

	#Use multiprocessing module to run parallel processes for each chromosome
	pool = mp.Pool()
	find_union_islands_partial = partial(find_union_islands,no_control,temp_dir_1,temp_dir_2)
	pool.map(find_union_islands_partial, chroms)
	pool.close()

	outfile_name = (args.treatment_file[0].replace('.bed','')+'-vs-'+args.treatment_file[1]+'-W'+str(args.window_size)+
					'-G'+str(args.gap_size)+'-E'+str(args.e_value)+'-union.island')
	outfile_path = args.output_directory+'/'+outfile_name

	with open(outfile_path,'w') as outfile:
		for chrom in chroms:
			union_island_list=np.load(chrom+'_union_output.npy')
			for island in union_island_list:
				output_line = island[0]+'\t'+str(island[1])+'\t'+str(island[2])+'\t'+str(island[3])+'\n'
				outfile.write(output_line)

if __name__ == "__main__":
	main(sys.argv)
