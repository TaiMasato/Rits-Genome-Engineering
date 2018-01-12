#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import itemgetter, attrgetter
from collections import defaultdict
import time
import re
import math
import string
import itertools
import unit1
import unit2



def AnnotatedUniquegRNADict(AnnotatedUniquegRNAFile):
	#key:AnnotatedUniquegRNA_list value:遺伝子名として辞書化
	annotated_unique_gRNA_list = [i for i in open(AnnotatedUniquegRNAFile,'r')]
	gene_name_unique_gRNA_list = [re.split('[=\t]',i.rstrip())[8] for i in open(AnnotatedUniquegRNAFile,'r')]
	annotated_unique_gRNA_dic = dict(zip(annotated_unique_gRNA_list ,gene_name_unique_gRNA_list ))
	return annotated_unique_gRNA_dic
	
	
def PickUpCentergRNA(annotated_unique_gRNA_dic , MapolyATG2stop_list):
	pass
	
	
	
		
def PickUpCentergRNA(add_center_gap_sort_AortedList_list):
	pregenename = ''
	#pregenename = add_center_gap_sort_AortedList_list[0][8].split('=')[0]
	gRNAs_per_gene_list = []
	tenchilist = []
	gRNAs_per_gene_list.append(add_center_gap_sort_AortedList_list[0])
	for x in xrange(len(add_center_gap_sort_AortedList_list)):
		[genename , centergap] = [add_center_gap_sort_AortedList_list[x][8].split('=')[0] , add_center_gap_sort_AortedList_list[x][9]]
		
		if pregenename == genename:
			gRNAs_per_gene_list.append(add_center_gap_sort_AortedList_list[x])
		if pregenename != genename:
			#====================１遺伝子で処理中
			
			#転地行列か
			tenchi = map(list,zip(*gRNAs_per_gene_list))[-1]
			#最小値のインデックス
			center_ind = tenchi.index(min(tenchi))
			list_len = len(gRNAs_per_gene_list)
			if list_len < 10:
				for num in xrange(list_len):
					f.write(str(gRNAs_per_gene_list[num][0]) + '\t' + str(gRNAs_per_gene_list[num][1]) + '\t' + str(gRNAs_per_gene_list[num][2]) + '\t' + str(gRNAs_per_gene_list[num][3]) + '\t' + str(gRNAs_per_gene_list[num][4]) + '\t' + str(gRNAs_per_gene_list[num][5]) + '\t' + str(gRNAs_per_gene_list[num][6]) + '\t' + str(gRNAs_per_gene_list[num][7]) + '\t' + str(gRNAs_per_gene_list[num][8]) + '\n' )
					#print gRNAs_per_gene_list[num][0],
					#print gRNAs_per_gene_list[num][1],
					#print gRNAs_per_gene_list[num][2],
					#print gRNAs_per_gene_list[num][3],
					#print gRNAs_per_gene_list[num][4],
					#print gRNAs_per_gene_list[num][5],
					#print gRNAs_per_gene_list[num][6],
					#print gRNAs_per_gene_list[num][7],
					#print gRNAs_per_gene_list[num][8]
			elif list_len >= 10:
				if (center_ind) < 5:
					for num in xrange(0 , 10):
						f.write(str(gRNAs_per_gene_list[num][0]) + '\t' + str(gRNAs_per_gene_list[num][1]) + '\t' + str(gRNAs_per_gene_list[num][2]) + '\t' + str(gRNAs_per_gene_list[num][3]) + '\t' + str(gRNAs_per_gene_list[num][4]) + '\t' + str(gRNAs_per_gene_list[num][5]) + '\t' + str(gRNAs_per_gene_list[num][6]) + '\t' + str(gRNAs_per_gene_list[num][7]) + '\t' + str(gRNAs_per_gene_list[num][8]) + '\n' )
						#print gRNAs_per_gene_list[num][0],
						#print gRNAs_per_gene_list[num][1],
						#print gRNAs_per_gene_list[num][2],
						#print gRNAs_per_gene_list[num][3],
						#print gRNAs_per_gene_list[num][4],
						#print gRNAs_per_gene_list[num][5],
						#print gRNAs_per_gene_list[num][6],
						#print gRNAs_per_gene_list[num][7],
						#print gRNAs_per_gene_list[num][8]
				elif (list_len - center_ind) < 5:
					for num in xrange(list_len - 10 , list_len):
						f.write(str(gRNAs_per_gene_list[num][0]) + '\t' + str(gRNAs_per_gene_list[num][1]) + '\t' + str(gRNAs_per_gene_list[num][2]) + '\t' + str(gRNAs_per_gene_list[num][3]) + '\t' + str(gRNAs_per_gene_list[num][4]) + '\t' + str(gRNAs_per_gene_list[num][5]) + '\t' + str(gRNAs_per_gene_list[num][6]) + '\t' + str(gRNAs_per_gene_list[num][7]) + '\t' + str(gRNAs_per_gene_list[num][8]) + '\n' )
						#print gRNAs_per_gene_list[num][0],
						#print gRNAs_per_gene_list[num][1],
						#print gRNAs_per_gene_list[num][2],
						#print gRNAs_per_gene_list[num][3],
						#print gRNAs_per_gene_list[num][4],
						#print gRNAs_per_gene_list[num][5],
						#print gRNAs_per_gene_list[num][6],
						#print gRNAs_per_gene_list[num][7],
						#print gRNAs_per_gene_list[num][8]
				elif (center_ind) >= 5 and (list_len - center_ind) >= 5:
					for num in xrange(center_ind - 5 , center_ind + 5):
						f.write(str(gRNAs_per_gene_list[num][0]) + '\t' + str(gRNAs_per_gene_list[num][1]) + '\t' + str(gRNAs_per_gene_list[num][2]) + '\t' + str(gRNAs_per_gene_list[num][3]) + '\t' + str(gRNAs_per_gene_list[num][4]) + '\t' + str(gRNAs_per_gene_list[num][5]) + '\t' + str(gRNAs_per_gene_list[num][6]) + '\t' + str(gRNAs_per_gene_list[num][7]) + '\t' + str(gRNAs_per_gene_list[num][8]) + '\n' )
						#print gRNAs_per_gene_list[num][0],
						#print gRNAs_per_gene_list[num][1],
						#print gRNAs_per_gene_list[num][2],
						#print gRNAs_per_gene_list[num][3],
						#print gRNAs_per_gene_list[num][4],
						#print gRNAs_per_gene_list[num][5],
						#print gRNAs_per_gene_list[num][6],
						#print gRNAs_per_gene_list[num][7],
						#print gRNAs_per_gene_list[num][8]			
	#====================１遺伝子で処理中
			gRNAs_per_gene_list = []
			gRNAs_per_gene_list.append(add_center_gap_sort_AortedList_list[x])
		pregenename = genename
		print time.time() - start,
		print 'sec'
		

if __name__ == "__main__":
	start = time.time()
	gff3File = '../assets/Mpolymorphav3.1.gene.gff3'
	AnnotatedUniquegRNAFile = '../assets/UniquegRNAlist_withAnnotation.txt'
	write_file = '../assets/PickUp_gRNAs_of_Center.txt'
	
	MapolyATG2stop_list = unit1.GetATG2stopPosition(gff3File)   #unit1
	MapolyATG2stop_list = unit2.MapolyATG2stopListDeleteSplicingVariant(MapolyATG2stop_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	ATG2stop_added_centerinfo_list = unit1.SarchCenterPositionandArea(MapolyATG2stop_list)   #unit1
	annotated_unique_gRNA_dic = AnnotatedUniquegRNADict(AnnotatedUniquegRNAFile) #辞書化
	PickUpCentergRNA(annotated_unique_gRNA_dic , ATG2stop_added_centerinfo_list)
	
	
	
	
	#sorted_annotated_unique_grna_list = unit2.AnnotatedgRNALsitDeleteSplicingVariant(sorted_annotated_unique_grna_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	#AddCenterGapAortedList(ATG2stop_added_centerinfo_list , sorted_annotated_unique_grna_list)
	
	#PickUpCentergRNA(add_center_gap_sort_AortedList_list)
	f = open(write_file , 'w')
	f.close
	
	print time.time() - start,
	print 'sec : Done'