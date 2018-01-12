#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
gffの情報を元に遺伝子の中央を調べるunit
PickUp_Center_of_gRNA17.10.04.pyで使う
"""

import time
import math
import string

def GetATG2stopPosition(gff3File):
	#gffからcdsのみを取り出し、各遺伝子ごとにATGからstopの位置をまとめて使いやすい形にしてreturn
	
	
	MapolyCDS_line = []
	MapolyCDS_list = []
	MapolyATG2stop_list = []
	start = 0
	end = float('inf')
	One_Before_gRNA_Name = 'genename'
	
	for gff_mapoly_line in open(gff3File, 'r'):
		if "##" in gff_mapoly_line:
			continue
		[reference, source, feature, start, end, score, strand, frame, attr] = gff_mapoly_line.split('\t')
		if feature == "CDS":
			attr = attr.split('.')
			attr = attr[0] + '.' + attr[1]
			MapolyCDS_line.extend([reference , attr[3:] , strand , start , end , ])
			MapolyCDS_list.append(MapolyCDS_line)   #['scaffold_998', 'Mapoly0998s0001.1', '+', '4625', '4861']
			MapolyCDS_line = []
	for x in xrange(len(MapolyCDS_list)):
		if MapolyCDS_list[x][1] != One_Before_gRNA_Name:
			MapolyATG2stop_list.append(MapolyCDS_list[x])
			start = 0
			end = float('inf')
		else:
			if MapolyCDS_list[x][2] == '+' and MapolyCDS_list[x][3] < start:
				MapolyATG2stop_list[-1][3] = MapolyCDS_list[x][3]
			elif MapolyCDS_list[x][2] == '-' and MapolyCDS_list[x][3] > start:
				MapolyATG2stop_list[-1][3] = MapolyCDS_list[x][3]
			if MapolyCDS_list[x][2] == '+' and end < MapolyCDS_list[x][4]:
				MapolyATG2stop_list[-1][4] = MapolyCDS_list[x][4]
			elif MapolyCDS_list[x][2] == '-' and end > MapolyCDS_list[x][4]:
				MapolyATG2stop_list[-1][4] = MapolyCDS_list[x][4]
				
		One_Before_gRNA_Name = MapolyCDS_list[x][1]
	MapolyATG2stop_list.sort(key=lambda x:x[1])
	return MapolyATG2stop_list
	

	
def SarchCenterPositionandArea(MapolyATG2stop_list):
	#MapolyATG2stop_listに配列の中央情報を加える。MapolyATG2stop_listの後ろに中央値、各遺伝子のうち中央として扱われる範囲のスタート位置、同じく終了位を新しい要素としてextend→return
	center_range = 40 #←%表記
	center_pos = ''
	start_center_range = ''
	end_center_range = ''
	ATG2stop_added_centerinfo_list = []

	for ATG2stop_line in MapolyATG2stop_list:
		[reference , attr , strand , start , end , ] = ATG2stop_line
		len_ATG2stop = (int(end) - int(start))
		center_pos = int(start) + ((len_ATG2stop)* 1/ 10)
		half_len_center_range = (float(len_ATG2stop) * float(center_range)/100)/2   #真ん中から片方向に対して20%の長さ
		start_center_range = center_pos - int(half_len_center_range)
		end_center_range = center_pos + int(half_len_center_range)
		ATG2stop_added_centerinfo_line =  ATG2stop_line + [center_pos , start_center_range , end_center_range]
		ATG2stop_added_centerinfo_list.append(ATG2stop_added_centerinfo_line)
	return ATG2stop_added_centerinfo_list
