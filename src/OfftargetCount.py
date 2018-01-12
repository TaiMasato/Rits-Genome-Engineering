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


def OTCounter():

    BowtieResultFile = "../assets/BowtieResult_with_MMposition_List.txt"

    OT_num = 0
    pre_gene = ""
    OTList_line = []
    OTList = []
    first = 0
    for BowtieResultFile_line in open(BowtieResultFile):
        BowtieResultFile_line = BowtieResultFile_line.strip()
        [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num] = BowtieResultFile_line.split('\t')
        if first == 0:
            pre_gene = GeneName
            first = 1
            continue
        if pre_gene != GeneName:
            OTList_line.extend([pre_gene,str(OT_num)])
            OTList.append(OTList_line)
            OTList_line = []
            pre_gene = GeneName
            OT_num = 0
        else:
            OT_num += 1
    OTList_line.extend([pre_gene,str(OT_num)])
    OTList.append(OTList_line)
    
    f = open("../assets/OTList.txt","w")
    csvWriter = csv.writer(f)
    """
    for x in xrange(len(OTList)):
        f.write(OTList(x))
    f.close()
    """

    for line in OTList:
        csvWriter.writerow(line)
    f.close()

if __name__ == '__main__':
    OTCounter()




print time.time() - start,
print 'sec'
