#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import collections

def check_3gRNAs(openfile):
	genecount = 0
	sufficient_gene_count = 0
	insufficient_gene_count = 0
	insufficient_gene_list = []
	gRNAlist = [re.split( '[;=\t]', i.rstrip()) for i in open(openfile , 'r')]
	tenchi_gRNAlist = []
	tenchi_gRNAlist = list(map(list,zip(*gRNAlist))) #転置行列
	count_dict = collections.Counter(tenchi_gRNAlist[8])
	for genename , Possession in count_dict.items():
		genecount +=1
		if Possession >= 3:
			sufficient_gene_count += 1
		elif Possession < 3:
			insufficient_gene_count += 1
			insufficient_gene_list.extend([genename])
	print genecount
	print '3本ある:',
	print sufficient_gene_count
	print '3本ない:',
	print insufficient_gene_count
	print '3本ある遺伝子数/遺伝子総数:',
	print (float(sufficient_gene_count)/(float(sufficient_gene_count)+float(insufficient_gene_count)))*100,
	print '%'
	return [sufficient_gene_count , insufficient_gene_count , insufficient_gene_list]
		
def check_insufficient_gene_length(insufficient_gene_list , gff3file):
	gfflist =  [re.split( '[=.\t]', i.rstrip()) for i in open(gff3file , 'r') if '#' not in i]
	print '足りてない遺伝子名 遺伝子長'
	
	
	for insufficient_gene_name in insufficient_gene_list:
		insufficient_gene_name = insufficient_gene_name.split('.')[0]
		for gff_line in gfflist:
			if 'gene' == gff_line[2]:
				if gff_line[11] == insufficient_gene_name:
					
					print gff_line[11],
					print int(gff_line[4]) - int(gff_line[3])
	
if __name__ == '__main__':
	openfile = '../assets/selected_gRNAs.txt'
	print openfile
	gff3file = '../assets/Mpolymorphav3.1.gene.gff3'
	sufficient_gene_count , insufficient_gene_count , insufficient_gene_list = check_3gRNAs(openfile)
	check_insufficient_gene_length(insufficient_gene_list , gff3file)