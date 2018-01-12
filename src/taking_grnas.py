#!/usr/bin/env python
# -*- coding: utf-8 -*-

from operator import itemgetter
import time
import re
import itertools
start = time.time()

def grna_sort_per_gene(openfile):
	sort_list = []
	for grna_line in open(openfile , 'r'):
		grna_line = grna_line.rstrip()
		[reference , sourse , feature , start , end , score , strand , frame , attr] = grna_line.split('\t')
		[grnaname , seq , mm , seed , nonseed , ot] = attr.split(';')
		ot = ot[10:]
		ot = "Offtarget="+str(ot.zfill(4))
		[genename , strand , position , grnano] = grnaname.split('=')
		sort_list.append([reference , sourse , feature , start , end , score , strand , frame , attr , genename , seq , mm , seed , nonseed , ot])
	sort_list.sort(key=itemgetter(14))
	sort_list.sort(key=itemgetter(11,12 ,13) , reverse = True)
	sort_list.sort(key=itemgetter(9))
	return sort_list

def pickup3grnas(sort_list,gRNA_taking_num):
	sort_list_withOT = []
	sort_list_withoutOT = []
	withoutOTlist_loop = 0
	gRNA_Taking_List = []
	pre_gene = sort_list[0][9]
	
	for x in xrange(len(sort_list)):
		[reference , sourse , feature , start , end , score , strand , frame , attr , genename , seq , mm , seed , nonseed , ot] = sort_list[x]
		if pre_gene == genename:
			if ot != "Offtarget=0000":
				sort_list_withOT.append(sort_list[x])
			elif ot == "Offtarget=0000":
				sort_list_withoutOT.append(sort_list[x])
		elif pre_gene != genename:
			
			withoutOTlist_loop = int(gRNA_taking_num)
			if len(sort_list_withoutOT) < int(gRNA_taking_num):  #w/oOTで足りるだけgRNAを取り出す
				withoutOTlist_loop = len(sort_list_withoutOT)
			for i in xrange(withoutOTlist_loop):
				[referencet_withoutOT , sourset_withoutOT , featuret_withoutOT , startt_withoutOT , endt_withoutOT , scoret_withoutOT , strandt_withoutOT , framet_withoutOT , attrt_withoutOT , genenamet_withoutOT , seqt_withoutOT , mmt_withoutOT , seedt_withoutOT , nonseedt_withoutOT , ott_withoutOT] = sort_list_withoutOT[i]

				f.write(str(referencet_withoutOT) + '\t' + str(sourset_withoutOT) + '\t' + str(featuret_withoutOT) + '\t' + str(startt_withoutOT) + '\t' + str(endt_withoutOT ) + '\t' + str(scoret_withoutOT) + '\t' + str(strandt_withoutOT) + '\t' + str(framet_withoutOT) + '\t' + str(attrt_withoutOT) + '\n')  #リストを出力

			if withoutOTlist_loop >= gRNA_taking_num:   #w/oOTで作れなかった残りのgRNAを錬金
				continue
			rest_taking_loop = int(gRNA_taking_num) - int(withoutOTlist_loop)
			if len(sort_list_withOT) < rest_taking_loop:  #そもそも一遺伝子あたりのgRNAの数が少なすぎた場合
				rest_taking_loop = len(sort_list_withOT)
			for j in xrange(rest_taking_loop):
				[reference_withOT , sourse_withOT , feature_withOT , start_withOT , end_withOT , score_withOT , strand_withOT , frame_withOT , attr_withOT , genename_withOT , seq_withOT , mm_withOT , seed_withOT , nonseed_withOT , ot_withOT] = sort_list_withOT[j]

				f.write(str(reference_withOT) + '\t' + str(sourse_withOT) + '\t' + str(feature_withOT) + '\t' + str(start_withOT) + '\t' + str(end_withOT ) + '\t' + str(score_withOT) + '\t' + str(strand_withOT) + '\t' + str(frame_withOT) + '\t' + str(attr_withOT) + '\n')  #リストを出力

			sort_list_withOT = []
			sort_list_withoutOT = []
			if ot != "Offtarget=0000":
				sort_list_withOT.append(sort_list[x])
			elif ot == "Offtarget=0000":
				sort_list_withoutOT.append(sort_list[x])
		pre_gene = genename
	#以下帳尻合わせde最後の遺伝子を処理する
	withoutOTlist_loop = int(gRNA_taking_num)
	if len(sort_list_withoutOT) < int(gRNA_taking_num):  #w/oOTで足りるだけgRNAを取り出す
		withoutOTlist_loop = len(sort_list_withoutOT)
	for i in xrange(withoutOTlist_loop):
		[referencet_withoutOT , sourset_withoutOT , featuret_withoutOT , startt_withoutOT , endt_withoutOT , scoret_withoutOT , strandt_withoutOT , framet_withoutOT , attrt_withoutOT , genenamet_withoutOT , seqt_withoutOT , mmt_withoutOT , seedt_withoutOT , nonseedt_withoutOT , ott_withoutOT] = sort_list_withoutOT[i]
		#print sort_list_withoutOT[i][9]
		f.write(str(referencet_withoutOT) + '\t' + str(sourset_withoutOT) + '\t' + str(featuret_withoutOT) + '\t' + str(startt_withoutOT) + '\t' + str(endt_withoutOT ) + '\t' + str(scoret_withoutOT) + '\t' + str(strandt_withoutOT) + '\t' + str(framet_withoutOT) + '\t' + str(attrt_withoutOT) + '\n')  #リストを出力
	if withoutOTlist_loop >= gRNA_taking_num:   #w/oOTで作れなかった残りのgRNAを錬金
		rest_taking_loop = int(gRNA_taking_num) - int(withoutOTlist_loop)
	if len(sort_list_withOT) < rest_taking_loop:  #そもそも一遺伝子あたりのgRNAの数が少なすぎた場合
		rest_taking_loop = len(sort_list_withOT)
	for j in xrange(rest_taking_loop):
		[reference_withOT , sourse_withOT , feature_withOT , start_withOT , end_withOT , score_withOT , strand_withOT , frame_withOT , attr_withOT , genename_withOT , seq_withOT , mm_withOT , seed_withOT , nonseed_withOT , ot_withOT] = sort_list_withOT[j]
		f.write(str(reference_withOT) + '\t' + str(sourse_withOT) + '\t' + str(feature_withOT) + '\t' + str(start_withOT) + '\t' + str(end_withOT ) + '\t' + str(score_withOT) + '\t' + str(strand_withOT) + '\t' + str(frame_withOT) + '\t' + str(attr_withOT) + '\n')  #リストを出力

	
if __name__ == '__main__':
	openfile = '../assets/selected_GC_gRNAs5.txt'
	outputfile = '../assets/selected_gRNAs_central.txt'
	
	print "How many gRNAs do you want?"
	gRNA_taking_num = raw_input(">>>")
	sort_list = grna_sort_per_gene(openfile)
	f = open(outputfile , 'w')
	pickup3grnas(sort_list , gRNA_taking_num)
	f.close
	print time.time() - start,
	print 'sec : Done'
