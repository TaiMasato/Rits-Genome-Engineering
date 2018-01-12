#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import itemgetter, attrgetter
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

	for opened_annotated_unique_grna_line in open(AnnotatedUniquegRNAFile , 'r'):
		opened_annotated_unique_grna_line = opened_annotated_unique_grna_line.rstrip()
		[reference , source , feature , start , end , score , strand , frame , attribute] = opened_annotated_unique_grna_line.split('\t')
		annotated_unique_grna_line.extend([reference , source , feature , start , end , score , strand , frame , attribute])
		annotated_unique_grna_list.append(annotated_unique_grna_line)
		annotated_unique_grna_line = []
		sorted_annotated_unique_grna_list = annotated_unique_grna_list
	sorted_annotated_unique_grna_list.sort(key=lambda x:x[4])
	sorted_annotated_unique_grna_list.sort(key=lambda x:x[8].split('=')[0]) #各遺伝子ごとにgRNAposition昇順
	return sorted_annotated_unique_grna_list
	



def gRNAget_eachgene(sorted_annotated_unique_grna_list , ATG2stop_added_centerinfo_list):
	center_position_gap = 0
	
	
	for x,y in itertools.product(xrange(len(ATG2stop_added_centerinfo_list)) , xrange(len(sorted_annotated_unique_grna_list))):
		[center_info_gene_name , center_info_grna_center_position] = [ATG2stop_added_centerinfo_list[x][1] , ATG2stop_added_centerinfo_list[x][5]]
		[divided_gRNAs_gene_name , divided_gRNAs_start_position] = [sorted_annotated_unique_grna_list[y][8].split('=')[0] , sorted_annotated_unique_grna_list[y][3]]
		
		if center_info_gene_name != divided_gRNAs_gene_name:
			continue
		#print center_info_gene_name,
		#print sorted_annotated_unique_grna_list[y][8]
			#center_position_gap = abs(int(center_info_grna_center_position) - int(divided_gRNAs_start_position))
			#print center_position_gap
		
		





if __name__ == "__main__":
	start = time.time()
	gff3File = '../assets/Mpolymorphav3.1.gene.gff3'
	AnnotatedUniquegRNAFile = '../assets/UniquegRNAlist_withAnnotation.txt'
	write_file = '../assets/PickUp_gRNAs_of_Center.txt'

	MapolyATG2stop_list = unit1.GetATG2stopPosition(gff3File)   #unit1
	MapolyATG2stop_list = unit2.MapolyATG2stopListDeleteSplicingVariant(MapolyATG2stop_list)   #unit2 スプライシングバリアントのところ。必要に応じて消してね
	ATG2stop_added_centerinfo_list = unit1.SarchCenterPositionandArea(MapolyATG2stop_list)   #unit1
	sorted_annotated_unique_grna_list = SortAnnotatedUniquegRNA(AnnotatedUniquegRNAFile)
	sorted_annotated_unique_grna_list = unit2.AnnotatedgRNALsitDeleteSplicingVariant(sorted_annotated_unique_grna_list)   #unit2 スプライシングバリアントのところ。必要に応じて消してね
	for i in sorted_annotated_unique_grna_list:
		print i
	f = open(write_file , 'w')
	center_gRNAs_list = gRNAget_eachgene(sorted_annotated_unique_grna_list , ATG2stop_added_centerinfo_list)
	f.close

print time.time() - start,
print 'sec : Done'
