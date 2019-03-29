import os
import subprocess
import time

#SICER Modes
SICER = True
RECOGNICER = False

#Set one of these True
TIME = True
MEMORY = False
CORRECTNESS = False

test_dir =  '/nv/vol190/zanglab/jy2ma/sicer2'
data_path = test_dir+'/data'
new_sicer_result_path = test_dir+'/sicer_result/new'
old_sicer_result_path = test_dir+'/sicer_result/old'
old_sicer_path = test_dir+'/SICER_1'


files = ['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed',
        '45394_treat_rep1.bed','35394_treat_rep1.bed','35400_treat_rep1.bed','45393_treat_rep1.bed',
        '35393_treat_rep1.bed','45383_treat_rep1.bed','35389_treat_rep1.bed','35456_treat_rep1.bed',
        '45415_treat_rep1.bed','35417_treat_rep1.bed','35420_treat_rep1.bed','45384_treat_rep1.bed',
        '35408_treat_rep1.bed','35451_treat_rep1.bed','45404_treat_rep1.bed','35424_treat_rep1.bed',
        '35445_treat_rep1.bed']

gm12878_ctrl_group = set(['45210_treat_rep1.bed','36626_treat_rep1.bed','36629_treat_rep1.bed','45206_treat_rep1.bed',
        '36615_treat_rep1.bed','36646_treat_rep1.bed','45213_treat_rep1.bed','36579_treat_rep1.bed',
        '36609_treat_rep1.bed','45216_treat_rep1.bed','36581_treat_rep1.bed','36598_treat_rep1.bed',
        '45220_treat_rep1.bed','36572_treat_rep1.bed','36582_treat_rep1.bed','45212_treat_rep1.bed'])


def test_time():
    new_runtime_dict = {}
    old_runtime_dict = {}

    for f in files[0:2]:
        if f in gm12878_ctrl_group:
            c = 'GSM733742_GM12878_input.bed'
        else:
            c = 'GSM733780_K562_input.bed'

        #Execution of SICER2.0
        start = time.time()
        treat = data_path+'/'+f
        control = data_path+'/'+c
        new_output_dir = new_sicer_result_path + '/'+ f.replace('.bed','')
        subprocess.call(['sicer', 'SICER', '-t', 'treat', '-c', 'control' '-s', 'hg38', '-o', 'new_output_dir'])
        end = time.time()
        runtime = end-start
        new_runtime_dict[f] = runtime


        #Execution of old SICER
        start = time.time()
        sicer_sh = os.path.join(old_sicer_path, "SICER.sh")
        old_output_dir = old_sicer_result_path + '/'+ f.replace('.bed','')
        subprocess.call([sicer_sh, data_path, f, c, old_output_dir, 'hg38', '1', '200', '150', '0.74', '600', '0.01'])
        end = time.time()
        runtime = end-start
        old_runtime_dict[f] = runtime

    log = open("timing_log.txt", 'w+')
    for f in files:
        log.write("Timing Figure for "+f.replace('.bed','')+':\nOld: '
                    + str(old_runtime_dict[f]) +' s\tNew: '+ str(new_runtime_dict[f]) + ' s\n')
    log.close()


def test_memory():
    # new_mem_dict = {}
    # old_mem_dict = {}
    #
    # for f in files[0:2]:
    #     if f in gm12878_ctrl_group:
    #         c = 'GSM733742_GM12878_input.bed'
    #     else:
    #         c = 'GSM733780_K562_input.bed'
    #
    #     #Execution of SICER2.0
    #     start = time.time()
    #     new_output_dir = new_sicer_result_path + '/'+ f.replace('.bed','')
    #     subprocess.call("valgrind --tool=massif --massif-out-file=%s --sicer SICER -t %s -c %s -s hg38 -o %s", f, c, new_output_dir)
    #     end = time.time()
    #     runtime = end-start
    #     new_runtime_dict[f] = runtime
    #
    #
    #     #Execution of old SICER
    #     start = time.time()
    #     sicer_sh = os.path.join(old_sicer_path, "SICER.sh")
    #     old_output_dir = old_sicer_result_path + '/'+ f.replace('.bed','')
    #     subprocess.call("%s %s %s %s %s hg38 1 200 150 0.74 600 0.01", sicer_sh, data_path, f, c, old_output_dir)
    #     end = time.time()
    #     runtime = end-start
    #     old_runtime_dict[f] = runtime
    #
    # log = open("timing_log.txt", 'w+')
    # for f in files:
    #     log.write("Timing Figure for "+f.replace('.bed','')+':\nOld: '
    #                 + old_runtime_dict[f] +' s\tNew: '+ new_runtime_dict[f] + ' s\n')
    # log.close()
    print("tbd")

def test_correctness():
    print("tbd")


if __name__ == "__main__":
    print("===========================Running Test===========================")
    print("SICER:", SICER)
    print("RECOGNICER:",RECOGNICER)
    print("TIME:",TIME)
    print("MEMORY:",MEMORY)
    print("CORRECTNESS:",CORRECTNESS)

    if TIME:
        test_time()
    if MEMORY:
        test_memory()
    if CORRECTNESS:
        test_correctness()
