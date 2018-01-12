#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
splicing variant を考慮するかどうかを処理するunit
PickUp_Center_of_gRNA17.10.04.pyで使う
"""
import re
import time
import math
import string

def MapolyATG2stopListDeleteSplicingVariant(MapolyATG2stop_list):
	#MapolyATG2stop_listにおけるsplicing variantの考慮をしないようにする
	pregene = 'genename'
	variant_num = '1'
	MapolyATG2stop_del_splicing_variant_list = []

	for MapolyATG2stop_line in MapolyATG2stop_list:
		[reference , genename , strand , start , end ,] = MapolyATG2stop_line
		if pregene == genename.split('.')[0] and variant_num != genename.split('.')[1]:
			continue
		MapolyATG2stop_del_splicing_variant_list.append([reference , genename , strand , start , end ,])
		pregene = genename.split('.')[0]
	return MapolyATG2stop_del_splicing_variant_list

def AnnotatedgRNALsitDeleteSplicingVariant(sorted_unique_grna_list):
	genename = 'genename'
	variant_num = '1'
	sorted_unique_grna_del_splicing_variant_list = []
	
	for sorted_unique_grna_line in sorted_unique_grna_list:
		[attribute] = re.split('[.=]' , sorted_unique_grna_line[8]),
		if genename == attribute[0] and variant_num != attribute[1]:
			continue
		sorted_unique_grna_del_splicing_variant_list.append(sorted_unique_grna_line)
		genename = attribute[0]
		variant_num = attribute[1]
			
	#for i in sorted_unique_grna_del_splicing_variant_list:
		#print i	
	return sorted_unique_grna_del_splicing_variant_list
