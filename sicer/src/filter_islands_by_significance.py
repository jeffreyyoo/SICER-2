import re, os, sys, shutil
from math import *
from string import *

import multiprocessing as mp
import numpy as np
from functools import partial

from sicer.lib import GenomeData

def filter_by_fdr_SICER (args,chrom):
	file_name = args.treatment_file.replace('.bed','')+'_'+chrom+'_island_summary.npy'
	cutoff = args.false_discovery_rate
	summary_graph = np.load(file_name)
	summary_bed = []

	for line in summary_graph:
		if(line[7] <= cutoff):
			bed_line = (line[0],line[1],line[2],line[3])
			summary_bed.append(bed_line)

	np_summary_bed = np.array(summary_bed, dtype=object)
	np.save(file_name,np_summary_bed)

def filter_by_fdr_SICER_df (args, columnindex,chrom):
	file_name = chrom+'_union_island_summary.npy'
	cutoff = args.false_discovery_rate_df
	summary_graph = np.load(file_name)
	summary_bed = []

	for line in summary_graph:
		if(line[columnindex] <= cutoff):
			bed_line = (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12])
			summary_bed.append(bed_line)

	np_summary_bed = np.array(summary_bed, dtype=object)
	np.save(file_name,np_summary_bed)

def main(args,columnindex):
	chroms = GenomeData.species_chroms[args.species];
	total_island_count = 0
	total_read_count = 0

	df_call = (type(args.treatment_file)==list)  #Determines if this function was called by SICER or SICER-DF

	pool = mp.Pool(processes = min(mp.cpu_count(),len(chroms)))
	filtered_output=[]
	if(df_call):
		filter_by_fdr_partial = partial(filter_by_fdr_SICER_df,args,columnindex)
		filtered_output = pool.map(filter_by_fdr_partial, chroms)
	else:
		filter_by_fdr_partial = partial(filter_by_fdr_SICER,args)
		filtered_output = pool.map(filter_by_fdr_partial, chroms)

	outfile_name = ''
	if(df_call):
		outfile_name = (args.treatment_file[0].replace('.bed','')+'-W'+str(args.window_size)+'-G'+str(args.gap_size)+
						'-increased-islands-summary-FDR'+str(args.false_discovery_rate_df))
	else:
		outfile_name = (args.treatment_file.replace('.bed','')+'-W'+str(args.window_size)+'-G'
					+str(args.gap_size)+'-FDR'+str(args.false_discovery_rate)+'-island.bed')

	outfile_path = os.path.join(args.output_directory,outfile_name)
	with open(outfile_path,'w') as outfile:
		for chrom in chroms:
			island_file_name=''
			if(df_call):
				island_file_name = chrom+'_union_island_summary.npy'
			else:
				island_file_name = args.treatment_file.replace('.bed','')+'_'+chrom+'_island_summary.npy'
			island_list = np.load(island_file_name)
			for island in island_list:
				output_line = ''
				for i in range(0,len(island)):
					output_line += str(island[i])+'\t'
				output_line+='\n'
				outfile.write(output_line)
				total_island_count+=1
				total_read_count+=island[3]

	print ("Given significance",  str(args.false_discovery_rate), ", there are", total_island_count, "significant islands");
	return total_read_count

if __name__ == "__main__":
	main(sys.argv)
