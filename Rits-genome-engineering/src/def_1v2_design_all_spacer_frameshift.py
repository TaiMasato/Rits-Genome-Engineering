
#Frameshift用gRNAの設計。各遺伝子のCDS領域に対して指定されたPAMの上流20塩基をgff形式&multiple fastaとして出力
import re
import time

one_before_gene_name = 'genename'
grna_no_per_gene = 0

def search_gRNA_5_prime_pam(annotation_reference, annotation, feature, cds_start, cds_end, score, strand, frame, gene_name, DNA_to_be_processed, pam, fw_pam, rv_pam, pam_len, PamPosition, f3, f4):
	FwSpacerStart = PamPosition + pam_len		#gRNAの開始点　(Fw鎖)
	RvSpacerStart = PamPosition - 20		#gRNAの開始点　(Rv鎖)
	if pam == fw_pam:
		if strand == '+':
			gRNA_strand = "+"
		elif strand == '-':
			gRNA_strand = "-"
		grna_start_position = int(cds_start) + int (FwSpacerStart)		#ゲノム中のgRNA開始点
		FwCutPosition = int(cds_start) + int(FwSpacerStart) + 3				#ゲノム編集により切断される位置
		gRNA = DNA_to_be_processed[FwSpacerStart : FwSpacerStart + 20]		#gRNA配列
		if len(gRNA) >= 20:		#CDSの配列末端でgRNAが20ntに満たない場合の除外
			f3.write(str(annotation_reference) + '\t' + str(annotation) + '\t' + str(feature) + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+19) + '\t' + '.' + '\t' + str(gRNA_strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(FwCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + ';' + str(gRNA) + '\n' )
			f4.write('>' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(FwCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + '\n' + str(gRNA) + '\n' )
	elif pam == rv_pam:
		if strand == '+':
			gRNA_strand = "-"
		elif strand == '-':
			gRNA_strand = "+"
		grna_start_position = int(cds_start) + int (RvSpacerStart)
		RvCutPosition = int(cds_start) + int(PamPosition) - 3
		gRNA = DNA_to_be_processed[RvSpacerStart : RvSpacerStart + 20]
		if len(gRNA) >= 20:
			gRNA = gRNA.translate(str.maketrans('ATGC', 'TACG'))[::-1]			#逆鎖で表示
			f3.write(str(annotation_reference) + '\t' + str(annotation) + '\t' + str(feature) + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+19) + '\t' + '.' + '\t' + str(gRNA_strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(RvCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + ';' + str(gRNA) + '\n' )
			f4.write('>' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(RvCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + '\n' + str(gRNA) + '\n' )
	return gRNA


def search_gRNA_3_prime_pam(annotation_reference, annotation, feature, cds_start, cds_end, score, strand, frame, gene_name, DNA_to_be_processed, pam, fw_pam, rv_pam, pam_len, PamPosition, f3, f4):
	FwSpacerStart = PamPosition - 20		#gRNAの開始点　(Fw鎖)
	RvSpacerStart = PamPosition + pam_len		#gRNAの開始点　(Rv鎖)
	if pam == fw_pam:
		if strand == '+':
			gRNA_strand = "+" #+
		elif strand == '-':
			gRNA_strand = "+" #+
		grna_start_position = int(cds_start) + int (FwSpacerStart)		#ゲノム中のgRNA開始点
		FwCutPosition = int(cds_start) + int(PamPosition) - 3				#ゲノム編集により切断される位置
		gRNA = DNA_to_be_processed[FwSpacerStart : FwSpacerStart + 20]		#gRNA配列]
		if len(gRNA) >= 20:		#CDSの配列末端でgRNAが20ntに満たない場合の除外
			f3.write(str(annotation_reference) + '\t' + str(annotation) + '\t' + 'gRNA' + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+19) + '\t' + '.' + '\t' + str(gRNA_strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(FwCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + ';' + str(gRNA) + '\n' )
			f4.write('>' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(FwCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + '\n' + str(gRNA) + '\n' )
	elif pam == rv_pam:

		if strand == '+':
			gRNA_strand = "-"#-
		elif strand == '-':
			gRNA_strand = "-"
		grna_start_position = int(cds_start) + int (RvSpacerStart)
		RvCutPosition = int(cds_start) + int(RvSpacerStart) + 3
		gRNA = DNA_to_be_processed[RvSpacerStart : RvSpacerStart + 20]
		if len(gRNA) >= 20:
			gRNA = gRNA.translate(str.maketrans('ATGC', 'TACG'))[::-1]			#逆鎖で表示
			f3.write(str(annotation_reference) + '\t' + str(annotation) + '\t' + 'gRNA' + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+19) + '\t' + '.' + '\t' + str(gRNA_strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(RvCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + ';' + str(gRNA) + '\n' )
			f4.write('>' + str(gene_name) + '=' + str(gRNA_strand) + '=' + str(RvCutPosition).zfill(9) + '=' + str(grna_no_per_gene).zfill(4) + '\n' + str(gRNA) + '\n' )
	return gRNA

def search_match_pam_pos(annotation_reference, annotation, feature, cds_start, cds_end, score, strand, frame, gene_name, DNA_to_be_processed, pam, fw_pam, rv_pam, pam_len, f3, f4):
	global one_before_gene_name
	global grna_no_per_gene
	match = re.compile(pam)
	i = 0
	while i >= 0:
		matching = match.search(DNA_to_be_processed, i)
		if matching:
			if one_before_gene_name != gene_name:
				grna_no_per_gene = 1
				one_before_gene_name = gene_name
			PamPosition = matching.start()
			if pam == 'TTTV' or pam == 'VAAA':		#Cpf1用
				gRNA = search_gRNA_5_prime_pam(annotation_reference, annotation, feature, cds_start, cds_end, score, strand, frame, gene_name, DNA_to_be_processed, pam, fw_pam, rv_pam, pam_len, PamPosition, f3, f4)
			else:
				gRNA = search_gRNA_3_prime_pam(annotation_reference, annotation, feature, cds_start, cds_end, score, strand, frame, gene_name, DNA_to_be_processed, pam, fw_pam, rv_pam, pam_len, PamPosition, f3, f4)
			if len(gRNA) >= 20:
				grna_no_per_gene += 1		#遺伝子ごとのgRNA番号カウンター
			i = matching.start() + 1
		else:
			break

def MakegRNA(fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len , f3 , f4):
	for i in range(len(fasta_name_list)):
		main_genome_reference = fasta_name_list[i]
		main_genome_seq = fasta_seq_list[i]
		for annotation_line in annotation_list:
			[annotation_reference , annotation , feature , start , end , score , strand , frame , gene_name] = annotation_line
			annotation_reference_match = annotation_reference
			if str(annotation_reference_match) != str(main_genome_reference):		#fastaのseq名とgffのreference名前のうち, 一致するものを処理
				continue
			DNA_to_be_processed = main_genome_seq[int(start)-1:int(end)]		#cdsパートの抜き出し
			#print (fw_pam , rv_pam)
			search_match_pam_pos(annotation_reference , annotation , feature , start , end , score , strand , frame , gene_name , DNA_to_be_processed , fw_pam , fw_pam , rv_pam , pam_len  , f3 , f4)		#Fw鎖として処理
			search_match_pam_pos(annotation_reference , annotation , feature , start , end , score , strand , frame , gene_name , DNA_to_be_processed , rv_pam , fw_pam , rv_pam , pam_len  , f3 , f4)		#Rv鎖として処理

def design_gRNA_for_frameshift_func(output_file, output_file_for_bowtie, fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len):
	TimeMeasurement = time.time()       #処理時間の計算開始
	f3 = open(output_file, 'w')
	f4 = open(output_file_for_bowtie, 'w')
	MakegRNA(fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len , f3 , f4)
	f3.close()
	print(output_file + ' was written.')
	f4.close()
	print(output_file_for_bowtie + ' was written.')

	print ('gRNA Design Done!')
	print ((time.time() - TimeMeasurement), 'sec')
