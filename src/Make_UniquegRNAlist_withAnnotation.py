#-*- coding: utf-8 -*-
#!/usr/bin/env python
#日本語
"""
ファイル名等に()を使うとバグる
CTCTCTやAGAGAGAGなどのgRNAをBowtieにかけるとofftarget多すぎて結果が返ってこない
"""
import string
import csv
import re
import time
import commands
import sys
import collections
start = time.time()


def Compress_BowtieResult():

    BowtieResultFile = "../assets/Joined_AllgRNAquery_with_BowtieAnnotation.txt"
    pre_gene = ""
    BowtieResultCompressed_line = []
    BowtieResultCompressed_List = []
    flag = 1
    f = open("../assets/UniquegRNAlist_withAnnotation.txt","w")
    for BowtieResultFile_line in open(BowtieResultFile):
        BowtieResultFile_line = BowtieResultFile_line.strip()
# 出力の形式(2017.10.2. はゼニゴケのフォーマットに合わせた。)　000017	scaffold_1	CDS	Mapoly0230s0001.1=0017	+	2287	CCTCGCACTTACCCTCAATCTCT	1	1	0	79

        [refernce, source, feature, start, end, score, strand, frame, attribute, MM_num, Seed_num, Nonseed_num, OT_num] = BowtieResultFile_line.split('\t')
        GeneName = attribute.split(';')[0]
        if int(OT_num) == 0:
#            BowtieResultCompressed_line.extend([GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num, OT_num])
            attribute = attribute + ";MissMatch=" + str(MM_num) + ";Seed=" + str(Seed_num) + ";Nonseed=" + str(Nonseed_num) + ";OffTarget=" +str(OT_num)
            BowtieResultCompressed_line_str = '\t'.join([refernce, source, feature, start.zfill(7), end.zfill(7), score, strand, frame, attribute])
#	    print BowtieResultCompressed_line_str
            f.write(BowtieResultCompressed_line_str+"\n")
#            pre_gene = GeneName
#            BowtieResultCompressed_List.append(BowtieResultFile_line)
            BowtieResultCompressed_line = []
        if int(OT_num) != 0 and flag != 0 and pre_gene != GeneName:
            pre_gene = GeneName
            flag = 0
            continue
        if int(OT_num) != 0 and flag == 0:
            #BowtieResultCompressed_line.extend([GeneName, strand, scaffold, position, sequence,mismatch,MM_num,Seed_num,Nonseed_num,OT_num])
            #for element in BowtieResultCompressed_line:
                #f.write(str(element)+'\t')
            attribute = attribute + ";MissMatch=" + str(MM_num) + ";Seed=" + str(Seed_num) + ";Nonseed=" + str(Nonseed_num) + ";OffTarget=" +str(OT_num)
            BowtieResultCompressed_line_str = '\t'.join([refernce, source, feature, start.zfill(7), end.zfill(7), score, strand, frame, attribute])
            f.write(BowtieResultCompressed_line_str+"\n")
	    #print BowtieResultFile_line
            #BowtieResultCompressed_List.append(BowtieResultFile_line)
            BowtieResultCompressed_line = []
            flag = 1

#    print "output started"
#    for line in BowtieResultCompressed_List:
#	line_str = '\t'.join(line)
#        f.write(line_str+"\n")
    f.close()

if __name__ == '__main__':
    Compress_BowtieResult()




print time.time() - start,
print 'sec'
