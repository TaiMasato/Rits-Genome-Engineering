import re
import time

def compress_bowtie_result_func(BowtieResultFile, output_unique_gRNA_list_with_annotation):
    TimeMeasurement = time.time()       #処理時間の計算開始
    pre_gene = ""
    BowtieResultCompressed_line = []
    BowtieResultCompressed_List = []
    flag = 1
    f = open(output_unique_gRNA_list_with_annotation,"w")
    for BowtieResultFile_line in open(BowtieResultFile):
        BowtieResultFile_line = BowtieResultFile_line.strip()

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
    f.close()

    print ('Format Conversion2 Done!')
    print ((time.time() - TimeMeasurement), 'sec')
