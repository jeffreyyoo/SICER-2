import os
import subprocess
import time

try:
    import compare
except:
    print("Error: Need \"compare.py\" module.")

#SICER Modes
SICER = True
RECOGNICER = False

#Set one of these True
COREvsTIME = True
core_count = 4
TIME = False
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

def test_core_time(core=25):
    new_runtime_dict = {}
    old_runtime_dict = {}

    log = open("timing_result.txt", 'w+')
    for f in files:
        if f in gm12878_ctrl_group:
            c = 'GSM733742_GM12878_input.bed'
        else:
            c = 'GSM733780_K562_input.bed'

        #Execution of SICER2.0
        treat = data_path+'/'+f
        control = data_path+'/'+c
        new_output_dir = new_sicer_result_path + '/'+ f.replace('.bed','')
        start = time.time()
        subprocess.call(['sicer', '-t', treat, '-c', control, '-s', 'hg38', '-o', new_output_dir, '--cpu', core,'--wig_output'])
        end = time.time()
        runtime = end-start
        new_runtime_dict[f] = runtime

        print((f1.replace('.bed','')+':\t'+str(new_runtime_dict[f]) + '\n'))
        log.write(f1.replace('.bed','')+':\t'+str(new_runtime_dict[f]) + '\n')


    log.close()


def test_correctness():
    new_runtime_dict = {}
    old_runtime_dict = {}
    output_str_lst =[]

    log = open("df_accuracy_result.txt", 'w+')
    time_log = open("df_timing_result.txt", 'w+')

    for i in range(0,len(files),2):
        f1 = files[i]
        f2 = files[i+1]
        if f1 in gm12878_ctrl_group:
            c1 = 'GSM733742_GM12878_input.bed'
        else:
            c1 = 'GSM733780_K562_input.bed'

        if f2 in gm12878_ctrl_group:
            c2 = 'GSM733742_GM12878_input.bed'
        else:
            c2 = 'GSM733780_K562_input.bed'

        #Execution of SICER2.0
        t1 = data_path+'/'+f1
        path_c1 = data_path+'/'+c1
        t2 = data_path+'/'+f2
        path_c2 = data_path+'/'+c2
        new_output_dir = new_sicer_result_path + '/'+ (f1.replace('.bed','')+'_and_'+f2.replace('.bed',''))
        start = time.time()
        subprocess.call(['sicer_df', '-t', t1, t2, '-c', path_c1, path_c2, '-s', 'hg38', '-o', new_output_dir, '--wig_output'])
        end = time.time()
        runtime = end-start
        new_runtime_dict[f1] = runtime

        #Execution of old SICER
        sicer_df_sh = os.path.join(old_sicer_path, "SICER-df.sh")
        old_output_dir = old_sicer_result_path + '/'+ (f1.replace('.bed','')+'-and-'+f2.replace('.bed',''))
        if not os.path.exists(old_output_dir):
            os.mkdir(old_output_dir)
        test_module_path = os.getcwd()
        os.chdir(old_sicer_result_path)
        start = time.time()
        subprocess.call(['sh',sicer_df_sh, f1, c1, f2, c2, '200', '600', '0.01', '0.01'])
        os.chdir(test_module_path)
        end = time.time()
        runtime = end-start
        old_runtime_dict[f1] = runtime

        result = compare.df_compare(new_output_dir, old_output_dir, f1,f2)
        output_str = f1.replace('.bed','')+'-'+f2.replace('.bed','')+' Result: '+str(result)
        output_str_lst.append(output_str)
        print(output_str)

        log.write(output_str+'\n')
        time_log.write("Timing Figure for "+f1.replace('.bed','')+':\nOld: '
                    + str(old_runtime_dict[f1]) +' s\tNew: '+ str(new_runtime_dict[f1]) + ' s\n')


    log.close()
    time_log.close()

def test_time():
    new_runtime_dict = {}
    old_runtime_dict = {}

    log = open("timing_result.txt", 'w+')
    for f in files:
        if f in gm12878_ctrl_group:
            c = 'GSM733742_GM12878_input.bed'
        else:
            c = 'GSM733780_K562_input.bed'

        #Execution of SICER2.0
        treat = data_path+'/'+f
        control = data_path+'/'+c
        new_output_dir = new_sicer_result_path + '/'+ f.replace('.bed','')
        start = time.time()
        subprocess.call(['sicer', '-t', treat, '-c', control, '-s', 'hg38', '-o', new_output_dir, '--wig_output'])
        end = time.time()
        runtime = end-start
        new_runtime_dict[f] = runtime


        #Execution of old SICER
        sicer_sh = os.path.join(old_sicer_path, "SICER.sh")
        old_output_dir = old_sicer_result_path + '/'+ f.replace('.bed','')
        if not os.path.exists(old_output_dir):
            os.mkdir(old_output_dir)
        test_module_path = os.getcwd()
        os.chdir(old_sicer_result_path)
        start = time.time()
        subprocess.call(['sh',sicer_sh, data_path, f, c, old_output_dir, 'hg38', '1', '200', '150', '0.74', '600', '0.01'])
        end = time.time()
        os.chdir(test_module_path)
        runtime = end-start
        old_runtime_dict[f] = runtime

        print(("Timing Figure for "+f.replace('.bed','')+':\nOld: '+ str(old_runtime_dict[f]) +' s\tNew: '+ str(new_runtime_dict[f]) + ' s\n'))

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


if __name__ == "__main__":
    print("===========================Running Test===========================")
    print("SICER:", SICER)
    print("RECOGNICER:",RECOGNICER)
    print("COREvsTIME:",COREvsTIME)
    print("TIME:",TIME)
    print("MEMORY:",MEMORY)
    print("CORRECTNESS:",CORRECTNESS)
    print("==================================================================")
    if TIME:
        test_time()
    if MEMORY:
        test_memory()
    if CORRECTNESS:
        test_correctness()
    if COREvsTIME:
        test_core_time(core=core_count)
