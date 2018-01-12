#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import itemgetter, attrgetter
from collections import defaultdict
import time
import math
import string
import itertools
import unit1
import unit2



def SortAnnotatedUniquegRNA(AnnotatedUniquegRNAFile):
	#AnnotatedUniquegRNAFileを遺伝子名ごとにgRNA開始位置の早い順に並び替えたものをreturn
	annotated_unique_grna_line = []
	annotated_unique_grna_list = []
	sorted_annotated_unique_grna_list = []
	
	annotated_unique_grna_list = [i.rstrip("\n\r").split('\t') for i in open(AnnotatedUniquegRNAFile , 'r')]
	#[reference , source , feature , start , end , score , strand , frame , attribute]
	sorted_annotated_unique_grna_list = annotated_unique_grna_list
	sorted_annotated_unique_grna_list.sort(key=lambda x:x[4])
	sorted_annotated_unique_grna_list.sort(key=lambda x:x[8].split('=')[0]) #各遺伝子ごとにgRNAposition昇順(itemgetterにしてre.split()にするよりこれの方が早い)
	return sorted_annotated_unique_grna_list

def AddCenterGapAortedList(ATG2stop_added_centerinfo_list , sorted_annotated_unique_grna_list):
	centerpos_gap = 0
	add_center_gap_sort_aortedList_list = []
	for x,y in itertools.product(xrange(len(ATG2stop_added_centerinfo_list)) , xrange(len(sorted_annotated_unique_grna_list))):
		[center_info_attr , center_info_grna_center_position] = [ATG2stop_added_centerinfo_list[x][1] , ATG2stop_added_centerinfo_list[x][5]]
		[divided_gRNAs_gene_name , divided_gRNAs_start_position] = [sorted_annotated_unique_grna_list[y][8].split('=')[0] , sorted_annotated_unique_grna_list[y][3]]
		if center_info_attr != divided_gRNAs_gene_name:
			continue
		center_position_gap = abs(int(center_info_grna_center_position) - int(divided_gRNAs_start_position))
		add_center_gap_sort_aortedList_list.append(sorted_annotated_unique_grna_list[y] + [str(center_position_gap).zfill(5)])
		print time.time() - start,
		print 'sec'
	return add_center_gap_sort_aortedList_list
		
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
	print time.time() - start,
	print 'sec'
	MapolyATG2stop_list = unit2.MapolyATG2stopListDeleteSplicingVariant(MapolyATG2stop_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	print time.time() - start,
	print 'sec'
	ATG2stop_added_centerinfo_list = unit1.SarchCenterPositionandArea(MapolyATG2stop_list)   #unit1
	print time.time() - start,
	print 'sec'
	sorted_annotated_unique_grna_list = SortAnnotatedUniquegRNA(AnnotatedUniquegRNAFile)
	sorted_annotated_unique_grna_list = unit2.AnnotatedgRNALsitDeleteSplicingVariant(sorted_annotated_unique_grna_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	print time.time() - start
	add_center_gap_sort_AortedList_list = AddCenterGapAortedList(ATG2stop_added_centerinfo_list , sorted_annotated_unique_grna_list)
	print time.time() - start
	
	PickUpCentergRNA(add_center_gap_sort_AortedList_list)
	
	
	#for i in :
		#print i
	
	
	
	
	f = open(write_file , 'w')
	f.close
	
	print time.time() - start,
	print 'sec : Done'