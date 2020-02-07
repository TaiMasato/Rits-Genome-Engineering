
import re
import time

def Join_Bowtiesummary_2query_func(UniquegRNAList, BowtieResultFile, output_all_gRNAquery_with_bowtie_annotation):
    TimeMeasurement = time.time()       #処理時間の計算開始
    f = open(output_all_gRNAquery_with_bowtie_annotation,"w")
    bowtie_result_dict = {}
    print ('# hit しなかった gRNA')
    for BowtieResultFile_line in open(BowtieResultFile, 'r'):
        BowtieResultFile_line = BowtieResultFile_line.rstrip()
        bowtie_result_dict.setdefault(BowtieResultFile_line.split('\t')[0], []).append(BowtieResultFile_line)
    for UniquegRNA_line in open(UniquegRNAList, 'r'):
        UniquegRNA_line = UniquegRNA_line.rstrip()
        if '#' in UniquegRNA_line:
            continue
        gRNA_name = UniquegRNA_line.split('\t')[8].split(';')[0]
        try:
            for bowtie_result_line in list(bowtie_result_dict[gRNA_name]):
                [GeneName, strand, scaffold, position, sequence, mismatch, MM_num, Seed_num, Nonseed_num, OT_num] = bowtie_result_line.rstrip().split('\t')
                f.write(str(UniquegRNA_line) + '\t' + str(MM_num) + '\t' + str(Seed_num) + '\t' + str(Nonseed_num) + '\t' + str(OT_num) + '\n')
        except KeyError:
            print (gRNA_name)
    f.close()

    print ('Format Conversion1 Done!')
    print ((time.time() - TimeMeasurement), 'sec')
