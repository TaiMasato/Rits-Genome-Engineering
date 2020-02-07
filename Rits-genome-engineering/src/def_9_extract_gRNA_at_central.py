import re
import time
def make_ORF_dict_per_gene(annotation_list):
	annotation_list.sort(key=lambda x:x[3])		#start位置でソート
	annotation_list.sort(key=lambda x:x[8])		#遺伝子名でソート
	annotation_dict_list = {}
	for annotation_line in annotation_list:		#複数要素を持つdictの作成
		if annotation_line[8] in annotation_line:
			annotation_dict_list.setdefault(annotation_line[8], []).append(tuple(annotation_line))		#key:遺伝子名, value: annotation line
	return annotation_dict_list

def make_center_list(annotation_dict_list, gene_id_list):
	ORF_center_list = []
	print ('KeyError list (gene_id_codeとannotation_dict間):')
	x=0
	for gene_id_code in gene_id_list:
		CDS_range_data = []
		try:
			for annotation_dict_line in list(annotation_dict_list[gene_id_code]):
				[reference, source, feature, start, end, score, strand, frame, gene_name] = list(annotation_dict_line)
				CDS_range_data.extend([int(start), int(end),])
		except KeyError:
			print ('KeyError:', gene_id_code)
			continue
		ORF_start = min(CDS_range_data)
		ORF_end = max(CDS_range_data)
		ORF_center = int(ORF_start) + int(abs(int(ORF_end)-int(ORF_start))/2)
		x+=1
		ORF_center_list.append([gene_name, ORF_center])
	return ORF_center_list

def make_gRNA_dict_list(unique_grna_with_annotated_file):		#def_10と同じコード　あとで共通化
	gRNA_dict_list = {}
	uniqe_grna_annotated_list = [j.rstrip().split() for j in open(unique_grna_with_annotated_file, 'r')]
	uniqe_grna_annotated_list.sort(key=lambda x:x[8].split('=')[2])		#cut_positionでソート
	uniqe_grna_annotated_list.sort(key=lambda x:x[8].split('=')[0])		#遺伝子名でソート
	for uniqe_grna_annotated_line in uniqe_grna_annotated_list:
		if str(uniqe_grna_annotated_line[8].split('=')[0]) in str(uniqe_grna_annotated_line):
			gRNA_dict_list.setdefault(uniqe_grna_annotated_line[8].split('=')[0], []).append(tuple(uniqe_grna_annotated_line))
	return gRNA_dict_list

def search_gRNA_center(ORF_center_list, gRNA_dict_list, write_file):
	pick_up_counter = 10		#中央付近から何本取り出すか？
	f = open(write_file, 'w')
	print ()
	print ('NOT hit gene names')
	for ORF_center_line in ORF_center_list:
		[gene_name, ORF_center] = ORF_center_line
		gRNA_with_center_dif_info_list = []
		try:
			for gRNA_per_gene in list(gRNA_dict_list[gene_name]):
				[reference, source, feature, start, end, score, strand, frame, attr] = list(gRNA_per_gene)
				cut_pos = list(gRNA_per_gene)[8].split('=')[2]
				gRNA_with_center_dif_info_list.append([reference, source, feature, start, end, score, strand, frame, attr, abs(int(ORF_center)-int(cut_pos))])
				#print (gene_name, abs(int(ORF_center)-int(cut_pos)))
			gRNA_with_center_dif_info_list.sort(key=lambda x:x[9])
			x = 0
			for pick_up_gRNA_line in gRNA_with_center_dif_info_list:
				[pick_up_reference, pick_up_source, pick_up_feature, pick_up_start, pick_up_end, pick_up_score, pick_up_strand, pick_up_frame, pick_up_attr, pick_up_center_diff] = pick_up_gRNA_line
				if x < pick_up_counter:
					#print (pick_up_attr,pick_up_center_diff, x)
					f.write(str(pick_up_reference) + '\t' + str(pick_up_source) + '\t' + str(pick_up_feature) + '\t' + str(pick_up_start) + '\t' + str(pick_up_end) + '\t' + str(pick_up_score) + '\t' + str(pick_up_strand) + '\t' + str(pick_up_frame) + '\t' + str(pick_up_attr) + '\n')
				x+=1
		except KeyError:
			print (gene_name)
		#print ()

def extract_gRNA_at_central_func(annotation_list, gene_id_list, unique_grna_with_annotated_file, write_file):
	TimeMeasurement = time.time()       #処理時間の計算開始
	annotation_dict_list = make_ORF_dict_per_gene(annotation_list)
	ORF_center_list = make_center_list(annotation_dict_list, gene_id_list)
	gRNA_dict_list = make_gRNA_dict_list(unique_grna_with_annotated_file)
	gRNA_center_info_list = search_gRNA_center(ORF_center_list, gRNA_dict_list, write_file)

	print ('Extract gRNA at Central Region Done!')
	print ((time.time() - TimeMeasurement), 'sec')
