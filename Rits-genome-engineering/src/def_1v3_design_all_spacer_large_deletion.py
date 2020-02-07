#Large deletion用gRNAの設計。各遺伝子の上下流領域に対して指定されたPAMの上流20塩基をgff形式&multiple fastaとして出力
import re
import time

def guide_RNA(gRNAstart, seq, pam, position, strand, rv_pam, pam_len, f3, f4): #position=各chromosomeの中での位置
	global designed_gRNAnum_upstream
	global designed_gRNAnum_downstream
	global gene_name
	global reference
	global stream
	global cut_gene

	#とりあえずNGANパターンを格納
	grna_start_position = position
	cut_position = position + 17 #Cas9によって切断される位置を格納
	seq = seq[gRNAstart:gRNAstart+20]

	if strand == "+":
		strand = "+"#+
	elif strand == "-":
		strand = "+"
	if pam == rv_pam:
		cut_position = position + 3
		seq = seq.translate(str.maketrans('ATGC','TACG'))[::-1]
		if strand == "+":
			strand = "-"
		elif strand == "-":
			strand = "-"
	if stream == 'u':
		f3.write(str(reference)+'\t'+'marker'+'\t'+'gRNA'+'\t'+str(grna_start_position)+'\t'+str(int(grna_start_position)+19)+'\t'+'.'+'\t'+str(strand)+'\t'+'.'+'\t'+str(gene_name)+"="+str(stream)+"="+str(strand)+"="+str(cut_position)+"="+str(designed_gRNAnum_upstream).zfill(4)+";"+str(seq)+str(cut_gene)+'\n')
		f4.write(">" + str(gene_name)+ "=" + str(stream) +"=" + str(strand)  + "=" + str(cut_position) + "=" + str(designed_gRNAnum_upstream).zfill(4) + "\n" + str(seq) + "\n")
		designed_gRNAnum_upstream += 1
	elif stream == 'd':
		f3.write(str(reference)+'\t'+'marker'+'\t'+'gRNA'+'\t'+str(grna_start_position)+'\t'+str(int(grna_start_position)+19)+'\t'+'.'+'\t'+str(strand)+'\t'+'.'+'\t'+str(gene_name)+"="+str(stream)+"="+str(strand)+"="+str(cut_position)+"="+str(designed_gRNAnum_downstream).zfill(4)+";"+str(seq)+str(cut_gene)+'\n')
		f4.write(">" + str(gene_name)+ "=" + str(stream) +"=" + str(strand)  + "=" + str(cut_position) + "=" + str(designed_gRNAnum_downstream).zfill(4) + "\n" + str(seq) + "\n")
		designed_gRNAnum_downstream += 1

def Select_For_or_Rev_gRNA(seq, pam, position, strand, fw_pam, rv_pam, pam_len, f3, f4):
	global gene_name
	global stream
	match = re.compile(pam)
	i = 0 #送られてきた範囲の塩基配列からNGAN(NTCN)の数を数える変数
	while i >= 0:
		matching = match.search(seq, i)

		if not matching:
			break
		global grnaList
		grnaList = []
		if pam == 'TTTC' or pam == 'GAAA':		#Cpf1用
			if pam == fw_pam: # PAM = NGG
				gRNAstart = (matching.start()) + int(pam_len)
				#if gRNAstart>=0:
				if gRNAstart <= len(seq)-20:
					guide_RNA(gRNAstart, seq, pam, position+gRNAstart, strand, rv_pam, pam_len, f3, f4)
			else:
				gRNAstart = (matching.start())-20
				#if gRNAstart <= len(seq)-20:
				if gRNAstart>=0:
					guide_RNA(gRNAstart, seq, pam, position+gRNAstart, strand, rv_pam, pam_len, f3, f4)
		else:
			if pam == fw_pam: # PAM = NGG
				gRNAstart = (matching.start())-20
				if gRNAstart>=0:
					guide_RNA(gRNAstart, seq, pam, position+gRNAstart, strand, rv_pam, pam_len, f3, f4)
					#print (stream,gene_name, gRNAstart, position, position+gRNAstart)
			else:
				gRNAstart = (matching.start()) + int(pam_len)
				if gRNAstart <= len(seq) - 20 - int(pam_len):
					guide_RNA(gRNAstart, seq, pam, position+gRNAstart, strand, rv_pam, pam_len, f3, f4)
					#print (stream,gene_name, gRNAstart, position, position+gRNAstart)
		i = matching.start() + 1

def digdown_gene(DNA_to_be_processed, genes_with_annotation_list, stream, strand, fw_pam, rv_pam, pam_len, f3, f4):
	global designed_gRNAnum_upstream
	global designed_gRNAnum_downstream
	global gene_name

	loop = 0#デフォルトで3回掘り下げる

	if stream == "u":
		pre_end = int(genes_with_annotation_list[3])-23 #遺伝子内を掘り下げていく時に一回前の掘り下げた3'末端の位置を格納する
		now_first = int(genes_with_annotation_list[3])+100 #遺伝子内を掘り下げていく時に今回掘り下げる5'末端の位置を格納する
		#print(stream, gene_name, pre_end)
		#print (genes_with_annotation_list)
		while designed_gRNAnum_upstream < 6:#ここまででgRNA候補の本数が6本未満の時
			Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , fw_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
			Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , rv_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
			pre_end = now_first
			now_first = now_first + 77 #77bpずつ下流にずらしていく
			loop += 1
			if loop == 3:#デフォルトで3回掘り下げる
				if designed_gRNAnum_upstream < 6:
				#ここで何かしらのflagを立てたい
					pass
				loop = 0
				break
	elif stream == "d":
		#print(stream, gene_name)

		pre_end = int(genes_with_annotation_list[4])-100
		now_first = int(genes_with_annotation_list[4])+23
		#print(stream, gene_name, pre_end)
		#print(stream, gene_name)
		#print (genes_with_annotation_list)
		while designed_gRNAnum_downstream < 6:#ここまででgRNA候補の本数が6本未満の時
			Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , fw_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
			Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , rv_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
			now_first = pre_end
			pre_end = pre_end - 77
			loop += 1
			if loop == 3:#デフォルトで3回掘り下げる
				if designed_gRNAnum_downstream < 6:
				#ここで何かしらのflagを立てたい
					pass
				loop = 0
				break
	"""
	loop = 0#デフォルトで3回掘り下げる
	if stream == "u":
		pre_end = int(genes_with_annotation_list[3])-23 #遺伝子内を掘り下げていく時に一回前の掘り下げた3'末端の位置を格納する
		now_first = int(genes_with_annotation_list[3])+100 #遺伝子内を掘り下げていく時に今回掘り下げる5'末端の位置を格納する
		#print(stream, gene_name)
		#print (genes_with_annotation_list)
	if stream == "d":
		pre_end = int(genes_with_annotation_list[4])-100
		now_first = int(genes_with_annotation_list[4])+23
		#print(stream, gene_name)
		#print (genes_with_annotation_list)
	while designed_gRNAnum < 6:#ここまででgRNA候補の本数が6本未満の時
		Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , fw_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
		Select_For_or_Rev_gRNA(DNA_to_be_processed[pre_end: now_first] , rv_pam, pre_end, strand, fw_pam, rv_pam, pam_len, f3, f4)
		if stream == "u":
			pre_end = now_first
			now_first = now_first + 77 #77bpずつ下流にずらしていく
		elif stream == "d":
			now_first = pre_end
			pre_end = pre_end - 77
		loop += 1
		if loop == 3:#デフォルトで3回掘り下げる
			if designed_gRNAnum < 6:
			#ここで何かしらのflagを立てたい
				pass
			loop = 0
			break
	"""


def chromosome_level_inquiry(genome_annotation_list, fasta_seq_list, fasta_name_list, fw_pam, rv_pam, pam_len, f3, f4):
	global designed_gRNAnum_upstream
	global designed_gRNAnum_downstream
	global each_gene_num
	global first
	global reference
	global chromosome_num
	global stream
	global gene_name
	global cut_gene

	designed_gRNAnum_upstream = 0 #各geneの5'側でそれぞれgRNAの本数を数える変数
	designed_gRNAnum_downstream = 0 #各geneの3'側でそれぞれgRNAの本数を数える変数
	each_gene_num = 0 #各chromosomeの遺伝子数を1つずつ足していく
	first = 0
	chromosome_num = 0 #chromosomeを1本ずつ足していく
	DNA_to_be_processed = ""
	cut_gene = ""
	pre_reference = ""

	for x in range(len(genome_annotation_list)):
		designed_gRNAnum_upstream = 0
		designed_gRNAnum_downstream = 0
		[reference , annotation , feature , start , end , score , strand , frame , gene_name] = genome_annotation_list[x]
		if reference != pre_reference:
		 	first = 0

		 	pre_end = ""
		 	pre_gene_name = ""
		 	pre_reference = reference
		 	for y in range(len(fasta_name_list)):
		 		if reference == fasta_name_list[y]:
		 			DNA_to_be_processed = fasta_seq_list[y]			#header1つ分のDNA配列
		stream = 'u'			#遺伝子の上流側
		temp = int(start) - 500			#遺伝子の開始点から500bp上流の間でgRNAを設計
		if temp < 0:
			temp = 0
		if first != 0:
			"""
			if str(gene_name)+str(stream) in not_enough_gene_list:
				cut_gene = str(" cutting ")+str(pre_gene_name)+str("d")
			else:
                cut_gene = ""
			"""
			if int(pre_end) > int(start)-500: #一つ前の遺伝子の終位置が次の遺伝子の開始位置-500bpより下流だったら
				#temp = int(end)
				temp = int(pre_end)		#前の遺伝子のendから今見ている遺伝子のstarの間でgRNAを設計
		Select_For_or_Rev_gRNA(DNA_to_be_processed[temp:int(start)], fw_pam, temp, strand, fw_pam, rv_pam, pam_len, f3, f4)
		Select_For_or_Rev_gRNA(DNA_to_be_processed[temp:int(start)], rv_pam, temp, strand, fw_pam, rv_pam, pam_len, f3, f4)
		if designed_gRNAnum_upstream < 6:		#設計できたgRNAが6本以下の時
			digdown_gene(DNA_to_be_processed, genome_annotation_list[x], stream, strand, fw_pam, rv_pam, pam_len, f3, f4)
			#pass
		pre_end = end
		pre_gene_name = gene_name
		pre_reference
		temp = int(end)+500		#遺伝子の終了点から500bp下流の間でgRNAを設計

		stream='d'			#遺伝子の下流側
		if x < len(genome_annotation_list)-1:		#リスト中で最後の遺伝子の時以外
			[next_reference , next_annotation , next_feature , next_start , next_end , next_score , next_strand , next_frame , next_gene_name] = genome_annotation_list[x+1]		#遺伝子の後ろ側に注目する
			if next_reference == reference:
				if int(next_start) < int(end)+500:
					temp = int(next_start)
		#print(stream, gene_name, pre_end)
		Select_For_or_Rev_gRNA(DNA_to_be_processed[int(end):temp], fw_pam, int(end), strand, fw_pam, rv_pam, pam_len, f3, f4)
		Select_For_or_Rev_gRNA(DNA_to_be_processed[int(end):temp], rv_pam, int(end), strand, fw_pam, rv_pam, pam_len, f3, f4)
		if designed_gRNAnum_downstream < 6:
			digdown_gene(DNA_to_be_processed, genome_annotation_list[x], stream, strand, fw_pam, rv_pam, pam_len, f3, f4)
			#pass

def design_gRNA_for_large_deletion_func(output_file, output_file_for_bowtie, fasta_name_list, fasta_seq_list, annotation_list, fw_pam , rv_pam , pam_len):
	TimeMeasurement = time.time()       #処理時間の計算開始
	f3 = open(output_file, "w")
	f3.write('# reference' + '\t' + 'sourse' + '\t' + 'feature' + '\t' + 'start' + '\t' +  'end' + '\t' + 'score' + '\t' + 'strand' + '\t' + 'frame' + '\t' + 'attribute' + '\n')
	f4 = open(output_file_for_bowtie, "w")
	chromosome_level_inquiry(annotation_list, fasta_seq_list, fasta_name_list, fw_pam, rv_pam, pam_len, f3, f4)
	f3.close()
	print(output_file + ' was written.')
	f4.close()
	print(output_file_for_bowtie + ' was written.')

	print ('gRNA Design Done!')
	print ((time.time() - TimeMeasurement), 'sec')
