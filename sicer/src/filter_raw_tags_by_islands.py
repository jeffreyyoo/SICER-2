#!/usr/bin/env python3

import re, os, sys, shutil
from math import *
from string import *
import bisect
import multiprocessing as mp
import numpy as np
from functools import partial

from sicer.lib import GenomeData

def tag_position(read, fragment_size):
    shift = int(round(fragment_size/2))
    if (read[5]=='+'):
        return int(read + shift)
    elif (read[5]=='-'):
        return read[2]-1-shift

def filter_tags_by_islands(file_name, fragment_size, chrom):
    island_list = np.load(file_name+'_'+chrom+'_island_summary.npy')
    read_list = np.load(file_name+'_'+chrom+'.npy')
    if(len(island_list)>0):
        island_start_list = []
        island_end_list = []
        for island in island_list:
            island_start_list.append(island[1])
            island_end_list.append(island[2])

        island_start_list.sort()
        island_end_list.sort()

        filtered_reads=[]
        for read in read_list:
            position = tag_position(read, fragment_size)
            if bisect.bisect_right(island_start_list, position) - bisect.bisect_left(island_end_list, position) == 1:
                filtered_reads.append(read)
        np_filtered_reads=np.array(filtered_reads,data=object)
        np.save(file_name+'_'+chrom+'.npy',np_filtered_reads)

def main(args):
    chroms = GenomeData.species_chroms[args.species];
    treatment_file = args.treatment_file.replace('.bed','')

    #Use multiprocessing to filter raw tags by islands in parallel processes
    pool = mp.Pool(processes = min(mp.cpu_count(),len(chroms)))
    filter_tags_by_islands_partial = partial(filter_tags_by_islands, treatment_file,args.fragment_size)
    pool.map(makeGraphFile_partial,chroms)
    pool.close()

    output_file_name = treatment_file+'-W'+str(args.window_size)
    if(args.subcommand == "SICER"):
        output_file_name +='-G'+str(args.gap_size)
    output_file_name += '-FDR'+args.false_discovery_rate+'-islandfiltered.bed'
    outfile_path = os.path.join(args.output_directory,output_file_name)
    with open(outfile_path,'w') as outfile:
        for chrom in chroms:
            filtered_bed = np.load(treatment_file+'_'+chrom+'.npy')
            for read in island_list:
                output_line = read[0]+'\t'+read[1]+'\t'+read[2]+'\t'+read[3]+'\t'+read[4]+'\t'+read[5]+'\n'
                outfile.write(output_line)

if __name__ == "__main__":
    main(sys.argv)