#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv
import string

def take(reference):
	for i in open('selected_gRNAs.txt', 'r'):
		i = i.rstrip()
		i = i.split('\t')
		
		if i[0] == reference:
			[grnaname , seq , mm , seed , nonseed, ot] = i[8].split(';')
			revseq = seq.translate(string.maketrans('ATGC','TACG'))[::-1]
			f.write(str(grnaname)+'-Fw' + '\t' + 'CTCG' +str(seq) + '\t' + '25nm' + '\t' + 'STD' + '\n' + str(grnaname)+ '-Rv' + '\t' + 'AAAC' +str(revseq) +'\t' + '25nm' + '\t' + 'STD' + '\n')		
		
if __name__ == '__main__':
	f = open('order.txt','w')
	take('scaffold_17')
	take('scaffold_18')
	take('scaffold_227')
	take('scaffold_230')
	take('scaffold_240')
	take('scaffold_250')
	take('scaffold_277')
	f.close()