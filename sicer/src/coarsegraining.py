#Author: Yiren Wang
#Edits made by Jeffrey Yoo

import re, os, sys, shutil
from math import *
from string import *
import bisect

import multiprocessing as mp
import numpy as np
from functools import partial

#From SICER package
from lib import GenomeData
from lib import Utility


'''version 8: 3-phase coarse graining, take the phase that has most 1 to next step. '''

def linreg(X, Y):
	"from Simple Recipes in Python http://www.phys.uu.nl/~haque/computing/WPark_recipes_in_python.html"
	if len(X) != len(Y):
		raise (ValueError, 'unequal length')
	N = len(X)
	Sx = Sy = Sxx = Syy = Sxy = 0.0
	for i in range(0,len(X)):
		x=X[i]
		y=Y[i]
		Sx = Sx + x
		Sy = Sy + y
		Sxx = Sxx + x*x
		Syy = Syy + y*y
		Sxy = Sxy + x*y
	det = Sxx * N - Sx * Sx
	if det != 0:
		return (Sxy * N - Sy * Sx)/det
	else:
		return 0
	#a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
	#return a


def is_list_sorted(List):
        """
        Check if sorted in ascending order.
        input is a list of pure numbers.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(List)-1):
                if List[index] > List[index + 1]:
                        sorted = 0;
        return sorted;


def start_list_correlation_r_rev(List, win, r, chrom_length):
	'''List must be sorted'''
	assert is_list_sorted(List) == 1
	x = List[0]%win
	d = int(r/win)
	SUMM = 0
	n = int((chrom_length - x)/win)
	if n - d > 0:
		a = [0] * n
		for item in List:
			i = int(item - x) // int(win)
			if i >= 0 and i < n:
				a[i] = 1
		for i in range(0, n - d):
			SUMM += a[i] * a[i + d]
		return float(SUMM)/float(n - d) - pow(float(sum(a))/float(len(a)),2)
	else:
		return 0.0


def start_list_correlation_function(List, win, chrom_length, name):
	xlist = []
	ylist = []
	#file = open("cr_"+name+"_"+str(win)+".txt", 'w')
	for i in range(0, min(3, int(chrom_length/win))):
		r = i * win
		c = start_list_correlation_r_rev(List, win, r, chrom_length)
		xlist.append(i)
		ylist.append(c)
		#file.write(str(i)+'\t'+str(c)+'\n')
	#file.close()
	return (xlist, ylist)


def correlation_length_fit(xlist, ylist):
	assert len(xlist) == len(ylist)
	loglist = []
	for i in range(0, len(ylist)):
		loglist.append(log(max(ylist[i], 0.000000000001)))
	a = linreg(xlist[1:],loglist[1:])
	if abs(a) > 0.000000000001:
		return -1.0/a
	else:
		return 1e12


def graining(List, win, step, score):
	'''
	1 step coarse graining, phase considered:
	List must be sorted!
	List (list) contains (start) coordinates of positive signals;
	win (int) is the window size in list, coarse graining will start from this resolution;
	step (int) is the number of windows in one graining unit;
	score (int) is the minimum number of positive elements in the graining unit to call the unit positive;
	output is a list of positive unit number in each graining step;
	'''
	result = []
	endlimit = List[-1]
	for p in range(0, step):
		tmp_result = []
		i = List[0] - p * win
		k = 0
		while i <= endlimit and k < len(List):
			j = i + step * win
			h = k
			while h <= (len(List)-1) and List[h] < j:
				h += 1
			n = h - k
			if n >= score:
				tmp_result.append(i)
			k = h
			i = j
		if len(tmp_result) > len(result):
			result = tmp_result
	return(result)


def coarsegraining(List, win_min, step, score, genome_length):
	if (is_list_sorted(List) != 1):
		List.sort()

	Length_list = []
	Length_list.append(len(List))  #number of eligible windows
	result_list = []
	result_list.append(List)	#list of start positions of eligible windows
	win = win_min
	while len(List) > 0:
		#(xlist, ylist) = start_list_correlation_function(List, win, genome_length)
		#print len(Length_list)-1, len(List)#, correlation_length_fit(xlist, ylist)
		List = graining(List, win, step, score)
		Length_list.append(len(List))
		if len(List) > 0:
			result_list.append(List)
		win = win * step
	return Length_list, result_list


def union_islands_to_list(islandlist, win):
	'''input islandlist and output list are both lists of BED island objects'''
	islandlist.sort(key=operator.attrgetter('start'));
	List = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end + 1 + win:
			List.append(current)
			current = compare
			i += 1
		else:
			current.end = max(current.end, compare.end)
			i += 1
	List.append(current)
	return List


def write_islandlist(List, win,chrom):
	'''input a start list and universal island width, output a islandlist of BED objects
	object.start = List[i]
	object.end = List[i] + win - 1'''
	output_list = []
	for start in List:
		output_list.append((chrom,start,int(start+win-1)))
	#output_list.sort(key=operator.attrgetter('start'))
	return output_list


def backstep(islandlist, List, win):
	'''one step trace back
		island[0]=chromosome
		island[1]=start position of island
		island[2]=end position of island'''
	#result_list = []
	#fine_islands = []
	addtional_islands = write_islandlist(List, win,chrom)
	for island in islandlist:
		start_left = (item[1] - win) in List
		start_right = item[1] in List
		if start_left and start_right:
			item[1] = item[1] - win
		elif (not start_left) and (not start_right):
			item[1] = item[1] + win
		end_left = (item[2] + 1 - win) in List
		end_right = (item[2] + 1) in List
		if end_left and end_right:
			item[2] = item[2] + win
		elif (not end_left) and (not end_right):
			item[2] = item[2] - win
		assert item[1] < item[2]
	return union_islands_to_list(islandlist + addtional_islands, win)


def traceback(List, win_min, step, level, genome_length, chrom):
	'''
	Input is a list of lists.
	'''
	win = win_min * pow(step, len(List)-1)
	islandlist = write_islandlist(List[-1], win,chrom)
	backlist = List[-1]
	(xlist, ylist) = start_list_correlation_function(backlist, win, genome_length, chrom)
	correlation_length = correlation_length_fit(xlist, ylist)
	#print len(backlist), correlation_length
	if len(List) > 1:
		(xlist, ylist) = start_list_correlation_function(List[-2], win/step, genome_length, chrom)
		correlation_length_next = correlation_length_fit(xlist, ylist)
		#print len(List[-2]), correlation_length_next
	i = 1
	while i < len(List)-level:
		backlist = List[-i-1]
		win = win/step
		if correlation_length > 1.0 and correlation_length_next >= correlation_length:
			break;
			#islands = islandlist
			#islandlist = backstep(islands, backlist, win)
			#if len(List) > i+1:
				#(xlist, ylist) = start_list_correlation_function(List[-i-2], win/step, genome_length, name)
				#print len(islandlist), correlation_length_fit(xlist, ylist)
		else:
			islandlist = write_islandlist(backlist, win,chrom)
			correlation_length = correlation_length_next
			if len(List) > i+1:
				(xlist, ylist) = start_list_correlation_function(List[-i-2], win/step, genome_length, chrom)
				correlation_length_next = correlation_length_fit(xlist, ylist)
				#print len(List[-i-2]), correlation_length_next
			else:
				correlation_length_next = 10000
		i += 1
	while i < len(List)-level:
		backlist = List[-i-1]
		#print len(islandlist)
		islands = islandlist
		islandlist = backstep(islands, backlist, win,chrom)
		win = win/step
		i += 1
	return islandlist

def filter_and_find_islands(args, min_tag_count, chrom):
	graph_file=args.treatment_file.replace('.bed','')+'_'+chrom+'_graph.npy'
	chrom_windows = np.load(graph_file)
	print_return = ''
	total_count_island = 0
	if(len(chrom_windows)>0):
		chrom_lengths = GenomeData.species_chrom_lengths[args.species][chrom]
		eligible_start_list = []

		for i in range(0,len(chrom_windows)):
			window = chrom_windows[i]
			read_count = window[3]
			if read_count >= min_tag_count:
				eligible_start_list.append(window[1])
		#print_return += "Coarse graining "+chrom+"\n"
		(result_list, island_list) = coarsegraining(eligible_start_list, args.window_size, args.step_size, args.step_score, chrom_lengths)
		#print_return += "Trace back "+chrom+"\n"
		islands = traceback(island_list, args.window_size, args.step_size, 0, chrom_lengths, chrom)

		if not(len(islands)>0):
			print_return +=chrom+" does not have any islands meeting the required significance"

		np_islands = np.array(islands,dtype=object)
		np.save(graph_file,np_islands)
		total_count_island=len(islands)

	return(total_count_island,print_return)

def main(args,read_count):

	print ("Coarse-graining approach to identify ChIP-Seq enriched domains:")
	print ("Species: ", args.species)
	print ("Window_size: ", args.window_size)
	print ("Coarse graining step: ", args.step_size)
	print ("Coarse graining score:", args.step_score)

	chroms = GenomeData.species_chroms[args.species]
	total_read_count = read_count
	print ("Total read count:", total_read_count)

	genome_length = sum (GenomeData.species_chrom_lengths[args.species].values());
	effective_genome_length = int(args.effective_genome_fraction * genome_length);

	average = float(total_read_count) * args.window_size/effective_genome_length;
	print ("Effective genome length: ", effective_genome_length)
	print ("window average:", average)

	min_tags_in_window = int(average) + 1
	print ("Minimum read count in a qualified window: ", min_tags_in_window)

	print ("Generate preprocessed data list")
	#read in the summary graph file
	#Use multiprocessing to find islands separately by chromosome
	pool = mp.Pool()
	filter_and_find_islands_partial = partial(filter_and_find_islands, args,min_tags_in_window)
	filtered_islands_result = pool.map(filter_and_find_islands_partial,chroms)
	pool.close()

	file_name = args.treatment_file.replace('.bed','')
	outfile_path = (args.output_directory+'/'+file_name+'-W'+str(args.window_size)
					+'-G'+str(args.gap_size)+'.cgisland')
	total_number_islands = 0
	path_to_filtered_graph = []
	with open(outfile_path,'w') as outfile:
		for i in range (0,len(filtered_islands_result)):
			filtered_chrom_graph = np.load(file_name+'_'+chroms[i]+'_graph.npy')
			total_number_islands+=filtered_islands_result[i][0]
			if(filtered_islands_result[i][1] != ""):
				print(filtered_islands_result[i][1])
			for island in filtered_chrom_graph:
				line = (island[0]+'\t'+str(island[1])+'\t'+str(island[2])
						+'\t'+'1\n')
				outfile.write(line)

	print ("Total number of islands: ", total_number_islands);



	# o = open(opt.out_file, 'w')
	# o.write('track type=bedGraph name=' + opt.out_file + '\n')
	# o.close()
	# SeparateByChrom.combineAllGraphFiles(chroms, ".islandstemp", opt.out_file)
	# SeparateByChrom.cleanup(chroms, ".islandstemp")
	# #else:
	# 	#print "input data error!"




if __name__ == "__main__":
	main(sys.argv)
