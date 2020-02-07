import re
import time
import datetime

def make_dict_reference_seq(FastaNameList , FastaSeqList):
	scaffold_namelist=[]
	scaffold_sequencelist=[]

	for scaffold_name in  open(FastaNameList,'r'): #fastaのreference名をリスト化
		if scaffold_name[0] == '>':
			scaffold_name = scaffold_name[1:]
		scaffold_name = scaffold_name.rstrip("\n\r")
		#print type(scaffold_name)
		scaffold_namelist.extend([scaffold_name])
	for scaffold_sequence in open(FastaSeqList,'r'): #fastaの配列をリスト化
		scaffold_sequence = scaffold_sequence.rstrip("\n\r")
		scaffold_sequencelist.extend([scaffold_sequence])
	dict_scaffold = dict(zip(scaffold_namelist , scaffold_sequencelist)) #辞書化
	return dict_scaffold

def pamcheck(dict_scaffold , BowtieResultFile , pam, f):
	BowtieResultList = []
	fw_pam = pam.translate(str.maketrans('N', '.'))
	rv_pam = pam.translate(str.maketrans('NATGC', '.TACG'))[::-1]
	Fw_PAM = re.compile(fw_pam)
	Rv_PAM = re.compile(rv_pam)
	for BowtieResultLine in open(BowtieResultFile , 'r'):
		BowtieResultLine = BowtieResultLine.rstrip("\n\r")
		BowtieResultLine = BowtieResultLine.split('\t')
		try:			#ミスマッチがなかった時の場合分け
			[gRNAname , strand , reference , position , seq , missmatch ] = BowtieResultLine
		except ValueError:
			[gRNAname , strand , reference , position , seq ] = BowtieResultLine
		if strand == '+':
			try:			#存在しないキーがある時の対処
				PAM_candidate = dict_scaffold[reference][int(position)+20:int(position)+20 +len(pam)]
#				print (gRNAname, PAM_candidate, missmatch)
			except KeyError:
				#print (reference, 'dose Not exist as key!' )
				continue
			match_fw_pam = Fw_PAM.search(PAM_candidate ,0)
			if match_fw_pam:
				#BowtieResultLine.append('True')
				for element in BowtieResultLine:
					f.write(str(element)+"\t")
				f.write("\n")   #+のPAM有り
		elif strand == '-':
			try:			#存在しないキーがある時の対処
				PAM_candidate = dict_scaffold[reference][int(position)-len(pam):int(position)]
#				print (gRNAname, PAM_candidate, missmatch)
			except KeyError:
				#print (reference, 'dose Not exist as key!' )
				continue
			match_rv_pam = Rv_PAM.search(PAM_candidate ,0)
			if match_rv_pam:
				#BowtieResultLine.append('True')
				for element in BowtieResultLine:
					f.write(str(element)+"\t")
				f.write("\n")  #-のPAM有り
		else:
			print ('Somthing wrong')

def pam_check_func(output_pam_check, output_file_fasta_name, output_file_fasta_seq, bowtie_result,inputted_pam):
	TimeMeasurement = time.time()       #処理時間の計算開始
	dict_scaffold = make_dict_reference_seq(output_file_fasta_name, output_file_fasta_seq)
	print (output_pam_check)
	f =open(output_pam_check,'w')
	pamcheck(dict_scaffold , bowtie_result ,inputted_pam, f)
	f.close()
	print ('Off-target Search Done!')
	print ((time.time() - TimeMeasurement), 'sec')
