#!/usr/bin/env python
#
# Authors: Chongzhi Zang, Weiqun Peng
#
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
from optparse import OptionParser
import operator
import bisect

from lib.GenomeData import *

plus = re.compile("\+");
minus = re.compile("\-");

def tag_position(sline, fragment_size):
	shift = int(round(fragment_size/2))
	if plus.match(sline[5]):
		return int(sline[1]) + shift
	elif minus.match(sline[5]):
		return int(sline[2]) - 1 - shift


def filter_tags_by_islands(chroms, islands, fragment_size):
	for chrom in chroms:
		if chrom in islands.keys():
			island_start_list = []
			island_end_list = []
			for item in islands[chrom]:
				island_start_list.append(item.start)
				island_end_list.append(item.end)
			island_start_list.sort()
			island_end_list.sort()
			bed_file = chrom + ".bed1"
			filtered_file = chrom +"_filtered.bed1"
			f = open(bed_file,'r')
			o = open(filtered_file, 'w')

			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					position = tag_position(sline, fragment_size)
					if bisect.bisect_right(island_start_list, position) - bisect.bisect_left(island_end_list, position) == 1:
						o.write('\t'.join(sline)+'\n')

			f.close()
			o.close()

def main(args,bed_file, island_file):

	if args.species in species_chroms.keys():
		chroms = species_chroms[args.species];
	else:
		print ("This species is not recognized, exiting");
		sys.exit(1);
	out_file = bed_file+"island_filtered.bed"
	islands = BED.BED(args.species, island_file, "BED3", 0);
	SeparateByChrom.separateByChrom(chroms, bed_file, '.bed1')
	filter_tags_by_islands(chroms, islands, args.fragment_size)
	final_output_file = out_file;
	final_output_file = SeparateByChrom.combineAllGraphFiles(chroms,
								   '_filtered.bed1',
								   final_output_file);
	SeparateByChrom.cleanup(chroms,'.bed1');
	SeparateByChrom.cleanup(chroms,'_filtered.bed1');


if __name__ == "__main__":
	main(sys.argv)
