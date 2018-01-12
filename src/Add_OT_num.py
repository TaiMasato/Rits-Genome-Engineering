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


def Add_OT_num():

    BowtieResultFile = open("../assets/BowtieResult_with_MMposition_List.txt")
    OTList = open("../assets/OTList.txt").readlines()
    OTList_line_num = 0
    BowtieResultAnnotatedList = []
    for BowtieResultFile_line in BowtieResultFile:
        BowtieResultFile_line = BowtieResultFile_line.strip()
        [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num] = BowtieResultFile_line.split('\t')
        for OTList_line_num in xrange(OTList_line_num,len(OTList)):
#            print "now calculating" + '\t'+ GeneName
            [OT_GeneName, OT_num] = OTList[OTList_line_num].split(',')
            if OT_GeneName == GeneName:
                BowtieResultFile_line = BowtieResultFile_line.strip() + "\t" + OT_num
                BowtieResultAnnotatedList.append(BowtieResultFile_line.strip())
                BowtieResultFile_line = []
                break
            if OT_GeneName != GeneName:
                OTList_line_num = OTList_line_num + 1 #OT_GeneNameとGeneNameが一致しない　⇒ そのOT_GeneNameに関しては、くまなくGeneNameを一致しおえた。
                continue




    f = open("../assets/BowtieResultAnnotatedFile_test.txt","w")
    for line in BowtieResultAnnotatedList:
        f.write(str(line)+"\n")
    f.close()

if __name__ == '__main__':
    Add_OT_num()



print "Add_OT_num.py has been done."
print time.time() - start,
print 'sec'
