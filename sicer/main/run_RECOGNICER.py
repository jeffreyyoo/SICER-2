# author:Jeffrey Yoo
# used python 3.5


import os
import shutil
import tempfile

curr_path = os.getcwd()

# From SICER Package
from sicer.src import remove_redundant_reads
from sicer.src import run_make_graph_file_by_chrom
from sicer.src import coarsegraining
from sicer.src import associate_tags_with_chip_and_control_w_fc_q
from sicer.src import filter_islands_by_significance
from sicer.src import make_normalized_wig
from sicer.src import filter_raw_tags_by_islands


def main(args, df_run=False):  # df_run indicates if run_RECOGNICER is being called by run_RECOGNICER_df function.

    # Checks if there is a control library
    control_lib_exists = True
    if (len(args.control_file) == 0):
        control_lib_exists = False

    # Ordinary version of SICER takes in one treatment file (and one control)
    if (type(args.treatment_file) == list):
        args.treatment_file = str(args.treatment_file[0])
    if (control_lib_exists and type(args.control_file) == list):
        args.control_file = str(args.control_file[0])

    try:
        dir_prefix = "SICER_" + str(os.getpid()) + "_"
        temp_dir = tempfile.mkdtemp(prefix=dir_prefix, dir=curr_path)
        # Change current working directory to temp_dir
        os.chdir(temp_dir)
    except:
        sys.exit(
            "Temporary directory required for SICER cannot be created. Check if directories can be created in %s." % args.input_directory)

    # Step 1: Remove redundancy reads in input file according to input threshold
    print("Preprocess the", args.treatment_file, "file to remove redundancy with threshold of",
          args.redundancy_threshold, "\n")
    treatment_file_path = os.path.join(args.input_directory, args.treatment_file)
    total_treatment_read_count = remove_redundant_reads.main(args, treatment_file_path)
    print('\n')

    # Step 2: Remove redundancy reads in control library according to input threshold
    total_control_reads = 0
    if (control_lib_exists):
        control_file_path = os.path.join(args.input_directory, args.control_file)
        print("Preprocess the", args.control_file, "file to remove redundancy with threshold of",
              args.redundancy_threshold, "\n")
        total_control_read_count = remove_redundant_reads.main(args, control_file_path)
        print('\n')

    # Step 3: Partition the genome in windows and generate graph files for each chromsome
    print("Partition the genome in windows and generate summary files \n")
    total_tag_in_windows = run_make_graph_file_by_chrom.main(args)
    print("\n")

    # Step4+5: Normalize and generate WIG file
    print("Normalizing graphs by total island filitered reads per million and generating summary WIG file \n")
    output_WIG_name = (args.treatment_file.replace('.bed', '') + "-W" + str(args.window_size) + "-normalized.wig")
    make_normalized_wig.main(args, output_WIG_name)

    # Step 6: Find condidate islands exhbiing clustering
    print("Finding candidate islands exhitiitng clustering \n")
    coarsegraining.main(args, total_tag_in_windows)
    print("\n")

    # Running SICER with a control library
    if (control_lib_exists):
        # Step 7
        print("Calculate significance of candidate islands using the control library \n")
        associate_tags_with_chip_and_control_w_fc_q.main(args, total_treatment_read_count, total_control_read_count)

        # Step 8:
        print("Identify significant islands using FDR criterion")
        significant_read_count = filter_islands_by_significance.main(args,
                                                                     7)  # 7 represents the ith column we want to filtered by
        print("Out of the ", total_treatment_read_count, " reads in ", args.treatment_file, ", ",
              significant_read_count, " reads are in significant islands")

    # Optional Outputs
    if (args.opt_output == 1):
        # Step 9: Filter treatment reads by the significant islands found from step 8
        print("Filter reads with identified significant islands...\n")
        filter_raw_tags_by_islands.main(args)

        # Step 10: Produce graph file based on the filtered reads from step 9
        print("Make summary graph with filtered reads...\n")
        run_make_graph_file_by_chrom.main(args)
        # Step 11: Produce Normalized WIG file
        print("Normalizing graphs by total island filitered reads per million and generating summary WIG file \n")
        output_WIG_name = (args.treatment_file.replace('.bed', '') + "-W" + str(args.window_size) + "-FDR" + str(
            args.false_discovery_rate) + "-islandfiltered-normalized.wig")
        make_normalized_wig.main(args, output_WIG_name)

    # Final Step
    if (df_run == True):
        return temp_dir
    else:
        print("Removing temporary directory and all files in it.")
        shutil.rmtree(temp_dir)
        print("End of SICER")
