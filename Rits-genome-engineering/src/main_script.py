import sys
import re
import time
import datetime

import message
import def_process_input_files
import def_1v2_design_all_spacer_frameshift
import def_1v3_design_all_spacer_large_deletion
import def_bowtie_process
import def_2_pam_check
import def_3_seed_check
import def_4_offtarget_count
import def_5_add_OT_num
import def_6_join_bowtie_summary_2query
import def_7_Make_Unique_gRNA_list_with_Annotation
import def_8_Restrict_gRNA
import def_9_extract_gRNA_at_central
import def_10_pick_out_gRNA



def main(argv, TimeMeasurement):
    try:
        [ desgin_pos, species_num,fasta, gff, gene_id, inputted_pam, restrict_gc, restrict_poly_t, extracted_central_region, sp_var, extracted_num] = argv     #frameshift用変数格納行
        print ('Species number [0; Arabidopsis thaliana/1; Marchantia polymorpha/2; Nicotiana benthamiana/3; Solanum lycopersicum/4; Populus_trichocarpa] : ', species_num)
        print ('gRNA Design Position [fs; Framshift/ld; Large Deletion] : ', desgin_pos)
        print ('Sequence File Name: ', fasta)
        print ('Annotation File Name: ', gff)
        print ('Gene ID File Name: ', gene_id)
        print ('PAM: ', inputted_pam)
        print ('Consider GC Content <y/n>: ', restrict_gc)
        print ('Consider Poly_T <y/n>: ', restrict_poly_t)
        print ('Extract gRNA from Central Region <y/n>: ', extracted_central_region)
        print (' Consider Splising Variant <y/n>: ', sp_var)
        print (' gRNA/gene <int/n>: ', extracted_num)
    except ValueError:
        try:
            [desgin_pos, species_num, fasta, gff, gene_id, inputted_pam, restrict_gc, restrict_poly_t, extracted_num] = argv     #frameshift用変数格納行
            print ('Species number [0; Arabidopsis thaliana/1; Marchantia polymorpha/2; Nicotiana benthamiana/3; Solanum lycopersicum/4; Populus_trichocarpa] : ', species_num)
            print ('gRNA Design Position [fs; Framshift/ld; Large Deletion] : ', desgin_pos)
            print ('Sequence File Name: ', fasta)
            print ('Annotation File Name: ', gff)
            print ('Gene ID File Name: ', gene_id)
            print ('PAM: ', inputted_pam)
            print ('Consider GC Content <y/n>: ', restrict_gc)
            print ('Consider Poly_T <y/n>: ', restrict_poly_t)
            print ('gRNA/gene <int/n>: ', extracted_num)
        except ValueError:
            message.messagefunc()       #input formatに適合していなければmessage表示
            return

    if int(species_num) == 0:
        species = 'Arabidopsis'
    elif int(species_num) == 1:
        species = 'Marchantia'
    elif int(species_num) == 2:
        species = 'N_benthamiana'
    elif int(species_num) == 3:
        species = 'S_lycopersicum'
    elif int(species_num) == 4:
        species = 'P_trichocarpa'

    if desgin_pos == 'fs':
        desgin = 'frameshift'

    elif desgin_pos == 'ld':
        desgin = 'large_deletion'

    if fasta.split('.')[-1] != 'fasta' and fasta.split('.')[-1] != 'fa' and fasta.split('.')[-1] != 'fas':    #input file1の拡張子がfasta形式かどうかの確認
        print ('File Extension Error')
        message.messagefunc()   #fasta形式でなければメッセージ表示
        return
    if gff.split('.')[-1] != 'gff3' and gff.split('.')[-1] != 'gff':    #input file2の拡張子がgff形式かどうかの確認
        print ('File Extension Error')
        message.messagefunc()   #gff形式でなければメッセージ表示
        return

    TimeMeasurement = time.time()       #処理時間の計算開始
    [year, month, day, hr, min, sec, msec] = re.split('[: \.-]',str(datetime.datetime.now()))      #日付の取得
    file_numbering = year[2:] + '_' + month + '_' + day + '_' + hr + '.' + min + '.' + sec      #year/month/day/hr:min:sec
    print (file_numbering)
    fasta = '../assets/input/' + fasta
    gff = '../assets/input/' + gff
    gene_id = '../assets/input/' + gene_id

    print ()
    print ('==== output file names ====')
    print ('fasta name list file: ', species + '_' + desgin + '_fasta_name_' + file_numbering + '.txt')
    output_file_fasta_name = '../assets/results/' + species + '_' + desgin + '_fasta_name_' + file_numbering + '.txt'
    print ('fasta seq list file: ', species + '_' + desgin + '_fasta_seq_' + file_numbering + '.txt')
    output_file_fasta_seq = '../assets/results/' + species + '_' + desgin + '_fasta_seq_' + file_numbering + '.txt'
    print ('output gff format file: ', species + '_' + desgin + '_allgRNAList_GFFformat_' + file_numbering + '.txt')
    output_file = '../assets/results/' + species + '_' + desgin + '_allgRNAList_GFFformat_' + file_numbering + '.txt'
    print ('output for bowite file: ', species + '_' + desgin + '_allgRNAList_for_bowtie_' + file_numbering + '.txt')
    output_file_for_bowtie = '../assets/results/' + species + '_' + desgin + '_allgRNAList_for_bowtie_' + file_numbering + '.txt'
    print ('bowtie result file: ', species + '_' + desgin + '_BowtieResultFile_' + file_numbering + '.txt')
    bowtie_result = '../assets/results/' + species + '_' + desgin + '_BowtieResultFile_' + file_numbering + '.txt'
    #print ('bowtie fixed file: ', species + '_' + desgin + '_BowtieResultFixedFile_' + file_numbering + '.txt')
    #output_file_bowtie_fixed_result = '../assets/results/' + species + '_' + desgin + '_BowtieResultFixedFile_' + file_numbering + '.txt'
    print ('pam checked file: ', species + '_' + desgin + '_BowtieResultPAMchecked' + file_numbering + '.txt')
    output_pam_check = '../assets/results/'+ species + '_' + desgin + '_BowtieResultPAMchecked_' + file_numbering + '.txt'
    print ('MMposition file: ', species + '_' + desgin + '_BowtieResult_with_MMposition_List_' + file_numbering + '.txt')
    output_MM_pos = '../assets/results/'+ species + '_' + desgin + '_BowtieResult_with_MMposition_List_' + file_numbering + '.txt'
    print ('Offtarget list file', species + '_' + desgin + '_OTList_' + file_numbering + '.txt')
    output_OTList = '../assets/results/'+ species + '_' + desgin + '_OTList_' + file_numbering + '.txt'
    print ('Bowtie result annotated file:', species + '_' + desgin + '_BowtieResultAnnotatedFile_' + file_numbering + '.txt')
    output_bowtie_result_annotated = '../assets/results/'+ species + '_' + desgin + '_BowtieResultAnnotatedFile_' + file_numbering + '.txt'
    print ('Joined all gRNAquery with bowtie annotation file: ',species + '_' + desgin + '_Joined_All_gRNAquery_with_Bowtie_Annotation_' + file_numbering + '.txt')
    output_all_gRNAquery_with_bowtie_annotation = '../assets/results/'+ species + '_' + desgin + '_Joined_All_gRNAquery_with_Bowtie_Annotation_' + file_numbering + '.txt'
    print ('unique gRNA list with annotation file: ', species + '_' + desgin + '_UniquegRNAlist_withAnnotation_'  + file_numbering +  '.txt')
    output_unique_gRNA_list_with_annotation = '../assets/results/'+ species + '_' + desgin + '_UniquegRNAlist_withAnnotation_'  + file_numbering +  '.txt'
    if restrict_gc == 'y' or restrict_poly_t =='y':
        print ('restricted gRNA list file: ', species + '_' + desgin + '_Restricted_gRNAs_'  + file_numbering +  '.txt')
        output_restrict_gRNA = '../assets/results/'+ species + '_' + desgin + '_Restricted_gRNAs_'  + file_numbering +  '.txt'
    if desgin_pos == 'fs' and extracted_central_region == 'y':
        print ('extracted gRNA at central CDS file: ', species + '_' + desgin + '_Extracted_gRNAs_at_Central_CDS_'  + file_numbering +  '.txt')
        output_extracted_gRNA_at_central_cds = '../assets/results/'+ species + '_' + desgin + '_Extracted_gRNAs_at_Central_CDS_'  + file_numbering +  '.txt'
    if extracted_num != 'n':
        print ('picked out gRNAs at central CDS file: ', species + '_' + desgin + '_picked_up_' + extracted_num + 'gRNAs_'  + file_numbering +  '.txt')
        output_picked_up_gRNAs_at_central_cds = '../assets/results/'+ species + '_' + desgin + '_picked_up_' + extracted_num + 'gRNAs_'  + file_numbering +  '.txt'



    print ()
    print ('==== Start processing ====')

    gene_id_list = def_process_input_files.get_gene_id_func(gene_id, desgin_pos, species_num,)     #種別操作あり
    if desgin_pos == 'fs':
        gene_id_list = def_process_input_files.process_splicing_variant_func(gene_id_list, sp_var)
    fasta_name_list, fasta_seq_list = def_process_input_files.fasta_file_func(fasta, output_file_fasta_name, output_file_fasta_seq)
    annotation_list = def_process_input_files.annotation_file_func(gff, desgin_pos, species_num, gene_id_list)     #種別操作あり
    fw_pam, rv_pam, pam_len = def_process_input_files.process_pam_func(inputted_pam)
#    def_process_input_files.mix_base_func(fw_pam, rv_pam)
    print ()
    print ('-- gRNA Design --')
         #gRNA設計、種別操作あり
    if desgin_pos == 'fs':
        def_1v2_design_all_spacer_frameshift.design_gRNA_for_frameshift_func(output_file, output_file_for_bowtie, fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len)   #frameshiftのgRNA設計
    elif desgin_pos == 'ld':
        def_1v3_design_all_spacer_large_deletion.design_gRNA_for_large_deletion_func(output_file, output_file_for_bowtie, fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len)   #largedeletionのgRNA設計


    ###bowtie処理まえにメモリ解放する必要ありそう。
    print ()
    print ()
    print ('-- Bowtie Process --')
    print ('# Bowtie log')
    def_bowtie_process.bowtie_process_func(species_num, output_file_for_bowtie, bowtie_result, fasta_name_list)     #bowtie内、種別操作あり
    print ()
    print ()
    print ('-- PAM Check --')
    def_2_pam_check.pam_check_func(output_pam_check, output_file_fasta_name, output_file_fasta_seq, bowtie_result,inputted_pam)
    print ()
    print ()
    print ('-- Miss Match Check --')
    def_3_seed_check.seed_check_func(output_pam_check, output_MM_pos)
    print ()
    print ()
    print ('-- Off-target Count --')
    def_4_offtarget_count.offtarget_count_func(output_MM_pos, output_OTList)
    print ()
    print ()
    print ('-- Add Off-target Number --')
    def_5_add_OT_num.add_OT_num_func(output_MM_pos, output_OTList, output_bowtie_result_annotated)
    print ()
    print ()
    print ('-- Format Conversion1 --')
    def_6_join_bowtie_summary_2query.Join_Bowtiesummary_2query_func(output_file, output_bowtie_result_annotated, output_all_gRNAquery_with_bowtie_annotation)
    print ()
    print ()
    print ('-- Format Conversion2 --')
    def_7_Make_Unique_gRNA_list_with_Annotation.compress_bowtie_result_func(output_all_gRNAquery_with_bowtie_annotation, output_unique_gRNA_list_with_annotation)

    used_file = output_unique_gRNA_list_with_annotation     #制限を行わないならこのファイルを次のプロセスに投げる
    if restrict_gc == 'y' or restrict_poly_t =='y':
        print ()
        print ()
        print ('-- Restricted gRNA List --')
        def_8_Restrict_gRNA.restrict_restrict_gRNA_func(restrict_gc, restrict_poly_t, output_unique_gRNA_list_with_annotation, output_restrict_gRNA)
        used_file = output_restrict_gRNA     #制限を行ったならこのファイルを次のプロセスに投げる

    if desgin_pos == 'fs' and extracted_central_region == 'y':
        print ()
        print ()
        print ('-- Extract gRNA at Central Region --')
        def_9_extract_gRNA_at_central.extract_gRNA_at_central_func(annotation_list, gene_id_list, used_file, output_extracted_gRNA_at_central_cds)
        used_file = output_extracted_gRNA_at_central_cds

    if extracted_num != 'n':
        print ()
        print ()
        print ('-- Extract <n>gRNAs/Gene --')
        def_10_pick_out_gRNA.pick_out_gRNA_func(desgin_pos, gene_id_list, species_num, used_file, output_picked_up_gRNAs_at_central_cds, extracted_num)

    print ()
    print ()

    print ('==== All process finished ====')
    print ((time.time() - TimeMeasurement), 'sec')

if __name__ == '__main__':
    TimeMeasurement = time.time()
    main(sys.argv[1:], TimeMeasurement)
