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


def Join_Bowtiesummary_2query():

    BowtieResultFile = open("../assets/BowtieResultAnnotatedFile_test.txt")
    UniquegRNAList = open("../assets/allgRNAList_GFFformat.txt").readlines()
    AllgRNA_num = 0
    BowtieResultAnnotatedList = []
    UniquegRNA_line = []
    f = open("../assets/Joined_AllgRNAquery_with_BowtieAnnotation.txt","w")
    for BowtieResultFile_line in BowtieResultFile:
        BowtieResultFile_line = BowtieResultFile_line.strip()
        [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num, OT_num] = BowtieResultFile_line.split('\t')
        for AllgRNA_num in xrange(AllgRNA_num,len(UniquegRNAList)):
            UniquegRNAList[AllgRNA_num].split('\t')
            UniquegRNA_line = UniquegRNAList[AllgRNA_num].split('\t')
            attribute = UniquegRNA_line[8].rstrip()
            UniquegRNA_line = UniquegRNA_line[:8] + [attribute]
            Unique_GeneName = attribute.split(';')[0]
            if Unique_GeneName == GeneName:
                BowtieResultFile_line = UniquegRNA_line + [MM_num] + [Seed_num] + [Nonseed_num] + [OT_num]
                BowtieResultFile_line_str = '\t'.join(BowtieResultFile_line)
                BowtieResultAnnotatedList.append(BowtieResultFile_line)
                f.write(BowtieResultFile_line_str+"\n")
                BowtieResultFile_line = []
                break
            if Unique_GeneName != GeneName:
                AllgRNA_num = AllgRNA_num + 1 #OT_GeneNameとGeneNameが一致しない　⇒ そのOT_GeneNameに関しては、くまなくGeneNameを一致しおえた。
                continue




#    f = open("../assets/BowtieResultCompressedFile_test.txt","w")
#    for line in BowtieResultAnnotatedList:
#        f.write(str(line)+"\n")
    f.close()

if __name__ == '__main__':
    Join_Bowtiesummary_2query()



print "Join_Bowtiesummary_2query.py has been done."
print time.time() - start,
print 'sec'
