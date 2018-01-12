#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
	
def SearchCentergRNAs(gRNAs_per_gene_list , ATG2stop_added_centerinfo_list):
	#2つのファイルを比較し遺伝子名が一致した時、(ATG2stop_added_centerinfo_listの中央値) - (gRNAs_per_gene_list)の値が宰相になったものが中央のgRNA。その両方向10本を拾い出力。
	pre_center_position_gap = float('inf')
	min_center_position_gap = 0
	center_gRNAs_list = []
	centerindex = -1
	pre_center_info_attr = ATG2stop_added_centerinfo_list[0][1]
	listlen = len(gRNAs_per_gene_list)
	
	for x , y in itertools.product(xrange(len(ATG2stop_added_centerinfo_list)) , xrange(len(gRNAs_per_gene_list))):
		[center_info_attr , center_info_grna_center_position] = [ATG2stop_added_centerinfo_list[x][1] , ATG2stop_added_centerinfo_list[x][5]]
		[divided_gRNAs_gene_name , divided_gRNAs_start_position] = [gRNAs_per_gene_list[y][8].split('=')[0] , gRNAs_per_gene_list[y][3]]
		print gRNAs_per_gene_list[y][8].split('=')
		if center_info_attr == divided_gRNAs_gene_name:
			centerindex += 1
			center_position_gap = abs(int(center_info_grna_center_position) - int(divided_gRNAs_start_position))
			if pre_center_position_gap < center_position_gap:
				break
			pre_center_position_gap = center_position_gap
	if len(gRNAs_per_gene_list) < 10:
		for num in xrange(len(gRNAs_per_gene_list)):
			#f.write(str(gRNAs_per_gene_list[y][0]) + '\t' + str(gRNAs_per_gene_list[y][1]) + '\t' + str(gRNAs_per_gene_list[y][2]) + '\t' + str(gRNAs_per_gene_list[y][3]) + '\t' + str(gRNAs_per_gene_list[y][4]) + '\t' + str(gRNAs_per_gene_list[y][5]) + '\t' + str(gRNAs_per_gene_list[y][6]) + '\t' + str(gRNAs_per_gene_list[y][7]) + '\t' + str(gRNAs_per_gene_list[y][8]) + '\n' )
			print 'O'
			#print gRNAs_per_gene_list[y][0],
			#print gRNAs_per_gene_list[y][1],
			#print gRNAs_per_gene_list[y][2],
			#print gRNAs_per_gene_list[y][3],
			#print gRNAs_per_gene_list[y][4],
			#print gRNAs_per_gene_list[y][5],
			#print gRNAs_per_gene_list[y][6],
			#print gRNAs_per_gene_list[y][7],
			#print gRNAs_per_gene_list[y][8]
			#print center_position_gap,
			#print centerindex,
			#print y
	elif len(gRNAs_per_gene_list) >= 10:
		if (centerindex) < 5:
			for num in xrange(-(centerindex) , 10 - (centerindex)):
				#f.write(str(gRNAs_per_gene_list[y-num][0]) + '\t' + str(gRNAs_per_gene_list[y-num][1]) + '\t' + str(gRNAs_per_gene_list[y-num][2]) + '\t' + str(gRNAs_per_gene_list[y-num][3]) + '\t' + str(gRNAs_per_gene_list[y-num][4]) + '\t' + str(gRNAs_per_gene_list[y-num][5]) + '\t' + str(gRNAs_per_gene_list[y-num][6]) + '\t' + str(gRNAs_per_gene_list[y-num][7]) + '\t' + str(gRNAs_per_gene_list[y-num][8]) + '\n' )
				print 'O'
				#print gRNAs_per_gene_list[y-num][0],
				#print gRNAs_per_gene_list[y-num][1],
				#print gRNAs_per_gene_list[y-num][2],
				#print gRNAs_per_gene_list[y-num][3],
				#print gRNAs_per_gene_list[y-num][4],
				#print gRNAs_per_gene_list[y-num][5],
				#print gRNAs_per_gene_list[y-num][6],
				#print gRNAs_per_gene_list[y-num][7],
				#print gRNAs_per_gene_list[y-num][8]
				#print y,
				#print y-num
		elif (listlen - centerindex) < 5:
			for num in xrange(-10 + (-10 + (listlen - (centerindex))) , listlen - (centerindex)):
				#f.write(str(gRNAs_per_gene_list[y-num][0]) + '\t' + str(gRNAs_per_gene_list[y-num][1]) + '\t' + str(gRNAs_per_gene_list[y-num][2]) + '\t' + str(gRNAs_per_gene_list[y-num][3]) + '\t' + str(gRNAs_per_gene_list[y-num][4]) + '\t' + str(RNAs_per_gene_list[y-num][5]) + '\t' + str(gRNAs_per_gene_list[y-num][6]) + '\t' + str(gRNAs_per_gene_list[y-num][7]) + '\t' + str(gRNAs_per_gene_list[y-num][8]) + '\n' )
				print 'O'
				#print gRNAs_per_gene_list[y-num][0]
				#print gRNAs_per_gene_list[y-num][0],
				#print gRNAs_per_gene_list[y-num][1],
				#print gRNAs_per_gene_list[y-num][2],
				#print gRNAs_per_gene_list[y-num][3],
				#print gRNAs_per_gene_list[y-num][4],
				#print gRNAs_per_gene_list[y-num][5],
				#print gRNAs_per_gene_list[y-num][6],
				#print gRNAs_per_gene_list[y-num][7],
				#print gRNAs_per_gene_list[y-num][8]
				#print y,
				#print y-num
		elif (centerindex) >= 5 and (listlen - centerindex) >= 5:
			for num in xrange(-5,5):
				#f.write(str(gRNAs_per_gene_list[y-num][0]) + '\t' + str(gRNAs_per_gene_list[y-num][1]) + '\t' + str(gRNAs_per_gene_list[y-num][2]) + '\t' + str(gRNAs_per_gene_list[y-num][3]) + '\t' + str(gRNAs_per_gene_list[y-num][4]) + '\t' + str(gRNAs_per_gene_list[y-num][5]) + '\t' + str(gRNAs_per_gene_list[y-num][6]) + '\t' + str(gRNAs_per_gene_list[y-num][7]) + '\t' + str(gRNAs_per_gene_list[y-num][8]) + '\n' )
				print 'O'
				#print gRNAs_per_gene_list[y-num][0],
				#print gRNAs_per_gene_list[y-num][1],
				#print gRNAs_per_gene_list[y-num][2],
				#print gRNAs_per_gene_list[y-num][3],
				#print gRNAs_per_gene_list[y-num][4],
				#print gRNAs_per_gene_list[y-num][5],
				#print gRNAs_per_gene_list[y-num][6],
				#print gRNAs_per_gene_list[y-num][7],
				#print gRNAs_per_gene_list[y-num][8]
				#print y,
				#print y-num

def gRNAget_eachgene(sorted_annotated_unique_grna_list , ATG2stop_added_centerinfo_list):
	#sorted_annotated_unique_grna_listファイルに記された遺伝子名を比較し、同じものをリスト化し１遺伝子ごとのgRNAリストを作成。SearchCentergRNAs()に１リストずつ飛ばす。
	gRNAs_per_gene_list = []
	pregenename = sorted_annotated_unique_grna_list[0][8].split('=')[0]
	
	for sorted_annotated_unique_grna_line in sorted_annotated_unique_grna_list:
		genename = sorted_annotated_unique_grna_line[8].split('=')[0]
		

		if pregenename == genename:
			gRNAs_per_gene_list.append(sorted_annotated_unique_grna_line)
		else:
			SearchCentergRNAs(gRNAs_per_gene_list , ATG2stop_added_centerinfo_list)
			gRNAs_per_gene_list = []
			gRNAs_per_gene_list.append(sorted_annotated_unique_grna_line)
		pregenename = genename
	SearchCentergRNAs(gRNAs_per_gene_list , ATG2stop_added_centerinfo_list)

if __name__ == "__main__":
	start = time.time()
	gff3File = '../assets/Mpolymorphav3.1.gene.gff3'
	AnnotatedUniquegRNAFile = '../assets/UniquegRNAlist_withAnnotation.txt'
	write_file = '../assets/PickUp_gRNAs_of_Center.txt'

	MapolyATG2stop_list = unit1.GetATG2stopPosition(gff3File)   #unit1
	MapolyATG2stop_list = unit2.MapolyATG2stopListDeleteSplicingVariant(MapolyATG2stop_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	ATG2stop_added_centerinfo_list = unit1.SarchCenterPositionandArea(MapolyATG2stop_list)   #unit1
	#for i in ATG2stop_added_centerinfo_list:
		#print i
	sorted_annotated_unique_grna_list = SortAnnotatedUniquegRNA(AnnotatedUniquegRNAFile)
	sorted_annotated_unique_grna_list = unit2.AnnotatedgRNALsitDeleteSplicingVariant(sorted_annotated_unique_grna_list)   #unit2 スプライシングバリアントのところ。考慮しない場合、必要に応じて消してね
	f = open(write_file , 'w')
	center_gRNAs_list = gRNAget_eachgene(sorted_annotated_unique_grna_list , ATG2stop_added_centerinfo_list)
	f.close

print time.time() - start,
print 'sec : Done'
