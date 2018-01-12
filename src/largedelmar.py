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
start = time.time()

GFF3='../assets/TAIR10_GFF3_genes_special_6_22.gff'
Genome_fasta = '../assets/TAIR10_chr_all.fas'

#GFF3 = '../assets/Mpolymorphav3.1.gene.gff3'
#Genome_fasta = '../assets/sca119.txt'

#----gfffileの読み込み
def get_gene_data_fromGFF(GFF3):
	
	genes_with_annotation_line = []
	genes_with_annotation_list = []
	chromosome_name_List = []
	prereference = 'referencename'
	genenum = 0# 各chromosomeのgene数
	geneNum_in_each_Chr = [] #各chromosomeのgene数を配列に格納する
	chromosome = 0 #chromosomeの数を数える変数
	
	for gff_line in open(GFF3):
		if "##" in gff_line: 
			continue
		gff_line = gff_line.rstrip()
		gff_line = gff_line.split('\t')
		[reference, source, feature, start, end, score, strand, frame, attr] = gff_line
		genename = re.split('[.;]',attr)[0][3:]
		
		if prereference == reference:
			if 'gene' == feature:
				genenum += 1
				genes_with_annotation_list.append([reference, source, feature, start, end, score, strand, frame, genename])
				
				#print reference,
				#print genename,
				#print feature,
				#print genenum
		elif prereference != reference:
			#print reference
			#print genenum
			chromosome_name_List.append(reference)			
			if 'gene' == feature:
				genenum = 1
				genes_with_annotation_list.append([reference, source, feature, start, end, score, strand, frame, genename])
				#print reference,
				#print genename,
				#print feature,
				#print genenum
		geneNum_in_each_Chr.append(genenum)
		prereference = reference
	print len(chromosome_name_List)

if __name__ == '__main__':
	
	#geneNum_in_each_Chr, chromosome, genes_with_annotation_list,chromosome_name_List = 
	get_gene_data_fromGFF(GFF3)
	
	print time.time() - start,
	print 'sec'