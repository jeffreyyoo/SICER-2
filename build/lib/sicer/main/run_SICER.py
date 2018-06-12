import tempfile
import os
import shutil

curr_path = os.getcwd() #Gets current working directory

#From SICER Package
from sicer.src import remove_redundant_reads
from sicer.src import run_make_graph_file_by_chrom
from sicer.src import find_islands_in_pr
from sicer.src import associate_tags_with_chip_and_control_w_fc_q
from sicer.src import filter_islands_by_significance
from sicer.src import make_normalized_wig

''' args: ArgumentParser object formed form command line parameters
	df_run: If df_run is true, then this instance of SICER is called by SICER-DF module.
			Default value is False.'''

def main(args, df_run = False):

	#Checks if there is a control library
	control_lib_exists = True
	if(args.control_file==""):
		control_lib_exists = False

	#Creates temporary directory to contain all intermediate files.
	dir_prefix = "SICER_"+str(os.getpid())+"_"
	temp_dir = tempfile.mkdtemp(prefix=dir_prefix,dir=curr_path)

	#Change current working directory to temp_dir
	os.chdir(temp_dir)

	#Step 1: Remove redundancy reads in input file according to input threshold
	#Output is the total number of reads retained. Represents size of library.
	print("Preprocess the", args.treatment_file, "file to remove redundancy with threshold of", args.redundancy_threshold,"\n")
	treatment_file_path = curr_path+'/'+args.treatment_file
	total_treatment_read_count= remove_redundant_reads.main(args, treatment_file_path)
	print('\n')

	#Step 2: Remove redundancy reads in control library according to input threshold
	total_control_reads=0
	if(control_lib_exists):
		control_file_path = curr_path+'/'+args.control_file
		print("Preprocess the", args.control_file, "file to remove redundancy with threshold of", args.redundancy_threshold,"\n")
		total_control_read_count = remove_redundant_reads.main(args,control_file_path)
		print('\n')

	#Step 3: Partition the genome in windows and generate graph files for each chromsome
	print("Partition the genome in windows and generate summary files \n")
	total_tag_in_windows = run_make_graph_file_by_chrom.main(args, args.treatment_file)
	print("\n")

	#Step4+5: Normalize and generate WIG file
	print("Normalizing graphs by total island filitered reads per million and generating summary WIG file \n")
	make_normalized_wig.main(args,args.treatment_file)

	#Step 6: Find condidate islands exhbiing clustering
	print("Finding candidate islands exhitiitng clustering \n")
	find_islands_in_pr.main(args,total_tag_in_windows)
	print("\n")

	#Running SICER with a control library
	if(control_lib_exists):
		#Step 7
		print("Calculate significance of candidate islands using the control library \n")
		significant_island_output = associate_tags_with_chip_and_control_w_fc_q.main(args,total_treatment_read_count,total_control_read_count)

		#Step 8:
		print("Identify significant islands using FDR criterion")
		significant_read_count = filter_islands_by_significance.main(args,7)  #7 represents the ith column we want to filtered by
		print("Out of the ",total_treatment_read_count," reads in ",args.treatment_file,", ",significant_read_count," reads are in significant islands")


	#Final Step
	if(df_run == False):
		print("Removing temporary directory and all files in it.")
		shutil.rmtree(temp_dir)
		print("End of SICER")
	else:
		return temp_dir
