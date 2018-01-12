#-*- coding: utf-8 -*-
#!/usr/bin/env python
#日本語

import string
import csv
import re
import time
import commands
import sys
import collections


#fastaname.txt
#fastaseq.txt
#BowtieResultFixedFile.txt

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

def pamcheck(dict_scaffold , BowtieResultFile , pam):
	BowtieResultList = []
	fw_pam = pam.translate(string.maketrans('N', '.'))
	rv_pam = pam.translate(string.maketrans('NATGC', '.TACG'))[::-1]
	Fw_PAM = re.compile(fw_pam)
	Rv_PAM = re.compile(rv_pam)
	for BowtieResultLine in open(BowtieResultFile , 'r'):
		BowtieResultLine = BowtieResultLine.rstrip("\n\r")
		BowtieResultLine = BowtieResultLine.split('\t')
		try:
			[gRNAname , strand , reference , position , seq , missmatch ] = BowtieResultLine
		except ValueError:
			[gRNAname , strand , reference , position , seq ] = BowtieResultLine
		if strand == '+':
			#print type(reference)
			PAM_candidate = dict_scaffold[reference][int(position)+20:int(position)+20 +len(pam)]
			match_fw_pam = Fw_PAM.search(PAM_candidate ,0)
			if match_fw_pam:
				#BowtieResultLine.append('True')
				for element in BowtieResultLine:
					f.write(str(element)+"\t")
				f.write("\n")   #+のPAM有り
		elif strand == '-':
			PAM_candidate = dict_scaffold[reference][int(position)-len(pam):int(position)]
			match_rv_pam = Rv_PAM.search(PAM_candidate ,0)
			if match_rv_pam:
				#BowtieResultLine.append('True')
				for element in BowtieResultLine:
					f.write(str(element)+"\t")
				f.write("\n")  #-のPAM有り
		else:
			print 'Somthing wrong'


if __name__ == '__main__':
	start = time.time()
	FastaNameList = '../assets/fastaname.txt'
	FastaSeqList = '../assets/fastaseq.txt'
	BowtieResultFile = '../assets/BowtieResultFixedFile.txt'
	output_file = '../assets/BowtieResultPAMchecked.txt'

	print 'What PAM?'
	print 'example : NGG'
	pam = raw_input('>>>')
	#pam='NGAN'
	dict_scaffold = make_dict_reference_seq(FastaNameList , FastaSeqList)
	f =open(output_file,'w')
	pamcheck(dict_scaffold , BowtieResultFile , pam)
	f.close()
	#dict_scaffold = {}
	#dict_scaffold = make_dict_reference_seq()
	print time.time() - start,
	print 'sec'
