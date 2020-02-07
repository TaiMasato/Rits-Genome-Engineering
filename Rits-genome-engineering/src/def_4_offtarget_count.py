import re
import csv
import time

def offtarget_count_func(output_MM_pos, output_OTList):
    TimeMeasurement = time.time()       #処理時間の計算開始
    OT_num = 0
    pre_gene = ""
    OTList_line = []
    OTList = []
    first = 0
    for BowtieResultFile_line in open(output_MM_pos):
        BowtieResultFile_line = BowtieResultFile_line.strip()
        [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num] = BowtieResultFile_line.split('\t')
        if first == 0:
            pre_gene = GeneName
            first = 1
            continue
        if pre_gene != GeneName:
            OTList_line.extend([pre_gene,str(OT_num).zfill(6)])
            OTList.append(OTList_line)
            OTList_line = []
            pre_gene = GeneName
            OT_num = 0
        else:
            OT_num += 1
    OTList_line.extend([pre_gene,str(OT_num).zfill(6)])
    OTList.append(OTList_line)
    f = open(output_OTList,"w")
    csvWriter = csv.writer(f)
    for line in OTList:
        csvWriter.writerow(line)
    f.close()

    print ('Off-target Count Done!')
    print ((time.time() - TimeMeasurement), 'sec')
