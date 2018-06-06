#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
import operator


def normalize_tag_count(input_file, ColumnIndex, total, output_file):
	"""
	Column index is 0 based 
	"""
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= (ColumnIndex + 1):
				sline[ColumnIndex] = str(float(sline[ColumnIndex])/total)
				outline = '\t'.join(sline) +'\n'
				outfile.write(outline)
	infile.close()
	outfile.close()


def total_counts(file, column):
	"""
	Column index is 0 based 
	"""
	total = 0.0
	infile = open(file,'r')
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= (column + 1):
				total += float(sline[column])
	infile.close()
	return total


def main(file, column, scaling_factor):
	total = total_counts(file, column) / scaling_factor
	normalize_tag_count(file, column, total, "normalized"+file)


if __name__ == "__main__":
	main(sys.argv)
