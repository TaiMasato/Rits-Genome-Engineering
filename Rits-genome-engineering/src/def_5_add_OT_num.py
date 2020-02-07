import re
import time

def add_OT_num_func(MM_pos, OTList, output_bowtie_result_annotated):
    TimeMeasurement = time.time()       #処理時間の計算開始
    OTList = open(OTList).readlines()
    OTList_line_num = 0
    BowtieResultAnnotatedList = []
    f = open(output_bowtie_result_annotated,'w')
    for BowtieResultFile_line in open(MM_pos, 'r'):
        BowtieResultFile_line = BowtieResultFile_line.strip()
        [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num] = BowtieResultFile_line.split('\t')
        for OTList_line_num in range(OTList_line_num,len(OTList)):
            [OT_GeneName, OT_num] = OTList[OTList_line_num].split(',')
            if OT_GeneName == GeneName:
                BowtieResultFile_line = BowtieResultFile_line.strip() + "\t" + OT_num
                BowtieResultAnnotatedList.append(BowtieResultFile_line.strip())
                f.write(str(BowtieResultFile_line.strip())+"\n")
                BowtieResultFile_line = []
                break
            if OT_GeneName != GeneName:
                OTList_line_num = OTList_line_num + 1 #OT_GeneNameとGeneNameが一致しない　⇒ そのOT_GeneNameに関しては、くまなくGeneNameを一致しおえた。
                continue
    f.close()

    print ('Add off-target number Done!')
    print ((time.time() - TimeMeasurement), 'sec')
