#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import re
import string
import collections
import csv
start = time.time()
BowtieResultFile = '../assets/BowtieResultPAMchecked.txt'

def Add_MMposition(BowtieResultFile):
	'''
	bowtieによって出力されたファイルをリスト化する。
	行と列を入れ替えたものを別のリストに格納し、各遺伝子に幾つのofftargetがあるかをカウント。
	offtargetの数を元のリストに新要素として格納、拡張する。
	'''
	BowtieResult_with_MMposition_line = []
	BowtieResult_ItemList = []
	BowtieResult_with_MMposition_List = []
	tr_BowtieResult_with_MMposition_List = []
	BowtieResult_with_MMcount_List = []
	seedarea = 12
	seed_num = 0
	nonseed_num = 0
	f = open("../assets/BowtieResult_with_MMposition_List.txt","w")
	for BowtieResultFile_line in open(BowtieResultFile, 'r'):   #ファイルの展開
		BowtieResultFile_line = BowtieResultFile_line.strip()
		BowtieResultFile_line = BowtieResultFile_line.split('\t')
		#[GeneName, strand, scaffold, position, sequence,mismatch]
		
		
		try:
			[GeneName, strand, scaffold, position, sequence,mismatch] = BowtieResultFile_line
		
		except ValueError:
			[GeneName, strand, scaffold, position, sequence] = BowtieResultFile_line
			mismatch =''
		
		sequence = sequence.strip()
		mismatch = mismatch.strip()
		MutationCounter = mismatch.count(':')

		'''
		Perfect matchのケース
		'''
		if MutationCounter == 0:
			BowtieResult_with_MMposition_line.extend([GeneName, strand, scaffold, position, sequence, mismatch,str(MutationCounter),str(seed_num),str(nonseed_num)])   #整頓したものを再リスト化(１次元)
			for element in BowtieResult_with_MMposition_line:
				f.write(str(element)+'\t')
			f.write(str('\n'))
			#BowtieResult_with_MMposition_List.append(BowtieResult_with_MMposition_line)
			 #再リスト化(２次元)
			BowtieResult_with_MMposition_line=[]
			continue
		'''
		Mutationが一つ以上のケース
		'''
		splitted_mismatch = re.split('[,:]',mismatch)
		for x in xrange(MutationCounter):
			if int(splitted_mismatch[2*x]) < seedarea:
				seed_num += 1
				continue
			nonseed_num += 1

		BowtieResult_with_MMposition_line.extend([GeneName, strand, scaffold, position, sequence,mismatch,str(MutationCounter),str(seed_num),str(nonseed_num)])   #整頓したものを再リスト化(１次元)
		for element in BowtieResult_with_MMposition_line:
			f.write(str(element)+'\t')
		f.write(str('\n'))
		#BowtieResult_with_MMposition_List.append(BowtieResult_with_MMposition_line)
		BowtieResult_with_MMposition_line=[]
		seed_num = 0
		nonseed_num = 0



	#csvWriter = csv.writer(f)
	"""
	for line in BowtieResult_with_MMposition_List:
		f.write(str(line))
	f.close()
	"""

if __name__ == '__main__':
	Add_MMposition(BowtieResultFile)   #bowtieの出力結果をリスト化すると共にMM情報を加える
	print time.time() - start,
	print 'sec'
