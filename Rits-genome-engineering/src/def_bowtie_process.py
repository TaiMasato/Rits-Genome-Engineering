#bowtie

import subprocess
#from memory_profiler import profile
import time

#@profile
def bowtie_check(species_num,output_file_for_bowtie, bowtie_result):
	TimeMeasurement = time.time()
	if species_num == '1':
		indexes = 'Mapoly'
	elif species_num == '0':
		indexes = 'TAIR10_chr_all'
	elif species_num == '2':
		indexes = 'Niben101'
	elif species_num == '3':
		indexes = 'ITAG4'
	elif species_num == '4':
		indexes = 'Ptri_v3'
	cmd='bowtie -a -v 3 ../assets/index/' + indexes + ' --suppress 6,7 -p 8 -f ' + output_file_for_bowtie +' >> ' + bowtie_result
	subprocess.call(cmd, shell=True)		#bowtieのコマンドを実行
	print ('Bowtie Done!')
	print ((time.time() - TimeMeasurement), 'sec')

def bowtie_result_fixer(bowtie_result, output_file_bowtie_fixed_result):
    f6 = open(output_file_bowtie_fixed_result,"w")
    with open(bowtie_result) as bowtie_f:
        for i, result_line in enumerate(bowtie_f):
            if result_line.startswith('#'):
                break
            result_line = str(result_line).rstrip()
            f6.write(result_line + '\n')
    f6.close()

def check_reference_pattern(fasta_name_list, output_file_bowtie_fixed_result):
    bowtie_result_reference = []
    with open(output_file_bowtie_fixed_result) as bowtie_fixed_f:
        for i, bowtie_fixed_line in enumerate(bowtie_fixed_f):
            bowtie_result_reference.extend([str(bowtie_fixed_line).rstrip().split('\t')[2]])
    bowtie_result_reference = list(set(bowtie_result_reference))
    print ('NOT Match Reference in Fasta Name')
    print (list(set(fasta_name_list)-set(bowtie_result_reference)))
    print ('NOT Match Reference in Bowtie Result')
    print (list(set(bowtie_result_reference)-set(fasta_name_list)))

#def bowtie_process_func(species_num, output_file_for_bowtie, bowtie_result, output_file_bowtie_fixed_result, fasta_name_list):
def bowtie_process_func(species_num, output_file_for_bowtie, bowtie_result, fasta_name_list):
    bowtie_check(species_num,output_file_for_bowtie, bowtie_result)
#    bowtie_result_fixer(bowtie_result, output_file_bowtie_fixed_result)
#    check_reference_pattern(fasta_name_list, output_file_bowtie_fixed_result)
