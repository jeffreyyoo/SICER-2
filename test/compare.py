import re
import sys
import math

control_file = "control.bed"
treatment_file = "treatment.bed"
output_files = ['treatment-W200-G600.scoreisland', 'treatment-W200-G600-FDR0.01-island.bed',
                'treatment-W200-G600-FDR0.01-islandfiltered.bed', 'treatment-W200-G600-FDR0.01-islandfiltered-normalized.wig',
                'treatment-W200-G600-islands-summary', 'treatment-W200-normalized.wig']


def check_islandsummary (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[5]= float(line1[5])
        line1[6]= float(line1[6])
        line1[7]= float(line1[7].replace('\n',''))
        line2 = re.split('\t',line2)
        line2[5]= float(line2[5])
        line2[6]= float(line2[6])
        line2[7]= float(line2[7].replace('\n',''))


        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[3]==line2[3]) 
                and (line1[4]==line2[4]) and (math.isclose(line1[5], line2[5], abs_tol=1e-7)) 
                and (math.isclose(line1[6], line2[6], abs_tol=1e-5)) and (math.isclose(line1[7], line2[7], abs_tol=1e-7)))

    return equal


def check_WIG (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        if(re.match("^track",line1)):
            equal = (line1 == line2)
        elif(re.match("^variableStep",line1)):
            equal = (line1 == line2)
        elif(line1!=""):
            line1 = re.split('\t',line1)
            line1[1]= float(line1[1].replace('\n',''))
            line2 = re.split('\t',line2)
            line2[1]= float(line2[1].replace('\n',''))

            equal == (line1[0]==line2[0]) and (math.isclose(line1[1], line2[1], abs_tol=0.01))

    return equal

def check_filteredbed(file1_name,file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line2 = re.split('\t',line2)

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[5]==line2[5]))

    return equal

def check_islandbed(file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[3]= line1[3].replace('\n','')
        line2 = re.split('\t',line2)
        line2[3]= line2[3].replace('\n','')

        equal = ((line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (line1[3]==line2[3]))

    return equal 


def check_scoreisland (file1_name, file2_name):
    file1 = open(file1_name,'r')
    file2 = open(file2_name,'r')
    equal = True
    for line1, line2 in zip(file1,file2):
        line1 = re.split('\t',line1)
        line1[3]= float(line1[3].replace('\n',''))
        line2 = re.split('\t',line2)
        line2[3]= float(line2[3].replace('\n',''))

        equal = (line1[0]==line2[0]) and (line1[1]==line2[1]) and (line1[2]==line2[2]) and (math.isclose(line1[3], line2[3], abs_tol=1e-6))

    return equal


def main():    
    chk_score_island = check_scoreisland(output_files[0], './expected_output/'+output_files[0])
    chk_island_bed = check_islandbed(output_files[1], './expected_output/'+output_files[1])
    chk_filtered_bed = check_filteredbed(output_files[2], './expected_output/'+output_files[2])
    chk_wig1 = check_WIG(output_files[3], './expected_output/'+output_files[3])
    chk_island_summary = check_islandsummary(output_files[4], './expected_output/'+output_files[4])
    chk_wig2 = check_WIG(output_files[5], './expected_output/'+output_files[5])

    return chk_score_island and chk_island_bed and chk_filtered_bed and chk_wig1 and chk_island_summary and chk_wig2

if __name__ == "__main__":
    result = main()
    if result:
        sys.exit(0)
    else:
        sys.exit(1)
