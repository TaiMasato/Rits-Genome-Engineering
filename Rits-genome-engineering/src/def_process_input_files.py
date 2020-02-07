#入力したファイルに対して
import re

def get_gene_id_func(gene_id_listfile, desgin_pos, species_num):
	gene_id_list = []
	for gene_id_listfile_line in open(gene_id_listfile, 'r'):
		if desgin_pos == 'fs':	#frameshift用gRNA設計の場合
			if species_num == '1':		#ゼニゴケのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}\.\d', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}\.\d', gene_id_listfile_line)	#ゼニゴケY染色体の正規表現
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d\.\d)', gene_id_listfile_line)
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '4':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
		elif desgin_pos == 'ld':		#large_deletion用gRNA設計の場合
			if species_num == '1':		#ゼニゴケのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}', gene_id_listfile_line)	#ゼニゴケY染色体の正規表現 #mayo
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d)', gene_id_listfile_line)
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '4':		#ポプラのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', gene_id_listfile_line)	#ゼニゴケ常染色体とX染色体の正規表現

		if gene_id:
			gene_id_name = gene_id.group()
			gene_id_list.extend([gene_id_name.strip()])
	gene_id_list = list(set(gene_id_list))
	gene_id_list.sort()
	return gene_id_list

def process_splicing_variant_func(gene_id_list, sp_var):
	if sp_var == 'n':			#gene_idからスプライシングバリアントの考慮をNoのときなくす
		gene_id_list_non_sp_list = []
		pre_gene_name = ''
		for gene_sp_var in gene_id_list:
			gene = gene_sp_var.split('.')[0]
			if str(pre_gene_name) != str(gene):
				gene_id_list_non_sp_list.extend([gene_sp_var])
			pre_gene_name = gene
		gene_id_list = []
		gene_id_list = gene_id_list_non_sp_list
	return gene_id_list

def fasta_file_func(genome_fasta, output_file_fasta_name , output_file_fasta_seq):
	fasta_name_list = []
	fasta_each_seq = ''
	fasta_seq_list = []
	f1 = open(output_file_fasta_name , 'w')  # 以下fastaの名前と配列リストの出力
	f2 = open(output_file_fasta_seq , 'w')

	opened_file = open(genome_fasta, 'r')		#ファイルopen
	for fasta_line in opened_file:		#一行ずつ取り出し
		fasta_line = fasta_line.rstrip()		#改行の削除
		if '>' in fasta_line:		#今forで見ている行がheaderの時
			if fasta_each_seq != '':		#条件文の構造上、最初にできるfasta_each_seqはから配列なので除去する
				fasta_seq_list.extend([fasta_each_seq.upper()])		#sORFごとに１行化した配列のリストに格納
				f2.write(str(fasta_each_seq) + '\n')
			header = fasta_line.split(' ')[0]
			fasta_name_list.extend([header[1:]])		#sORF_idをリストに格納(header除く)
			f1.write(str(header[1:]) + '\n')
			fasta_each_seq = ''		#次のheaderの配列が来るので古いheaderの配列を消して初期化
		else:		#今forで見ている行がseqの時
			fasta_each_seq += fasta_line		#50塩基ずつに区切られた塩基配列をsORFごとに合体する
	fasta_each_seq += fasta_line		#条件文の構造上、sORF_each_seqの最後はforで取り出せないのでここで超尻を合わせる
	fasta_seq_list.extend([fasta_each_seq.upper()])	#sORFごとに１行化した配列をリストに格納
	f2.write(str(fasta_each_seq).upper() + '\n')
	del fasta_each_seq 	#メモリ解放
	opened_file.close()		#ファイル閉じる
	f1.close
	f2.close
	return [fasta_name_list, fasta_seq_list]

def annotation_file_func(genome_annotation, desgin_pos, species_num, gene_id_list):
	annotation_dict_list = {}
	gene_name_for_search_list = []
	unique_gene_name_for_search_list = []
	for annotation_line in open(genome_annotation, 'r'):
		if "##" in annotation_line:
			continue
		[reference, source, feature, start, end, score, strand, frame, attr] = annotation_line.rstrip().split('\t')
		divided_attr = re.split('[.:;,-]',attr)
		if desgin_pos == 'fs':	#frameshift用gRNA設計の場合
			if species_num == '1':		#ゼニゴケのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}\.\d', attr)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}\.\d', attr)	#ゼニゴケY染色体の正規表現
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d\.\d)', attr)
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', attr)	#ゼニゴケ常染色体とX染色体の正規表
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', attr)	#ゼニゴケ常染色体とX染色体の正規表
			elif species_num == '4':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', attr)	#ゼニゴケ常染色体とX染色体の正規表
		elif desgin_pos == 'ld':		#large_deletion用gRNA設計の場合
			if species_num == '1':		#ゼニゴケのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}', attr)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}', attr)	#ゼニゴケY染色体の正規表現
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d)', attr)
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', attr)	#タバコ遺伝子の正規表現 #mayo
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', attr)	#トマト遺伝子の正規表現 #mayo
			elif species_num == '4':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', attr)	#ゼニゴケ常染色体とX染色体の正規表現 #mayo


		if gene_id:
			gene_idcode = gene_id.group()
		else:
			gene_idcode = "NONE"
		gene_name_for_search = gene_idcode
		annotation_gene_name = gene_idcode
		if species_num == '0':																#ナズナ帳尻合わせゾーン開始
			if reference != 'ChrC' and reference != 'ChrM':
				reference = reference[3:]
			elif reference == 'ChrM':
				reference = 'mitochondria'
			elif reference == 'ChrC':
				reference = 'chloroplast'													#ナズナ帳尻合わせゾーン終了
		line = [reference, source, feature, start, end, score, strand, frame, annotation_gene_name]
		gene_name_for_search_list.extend([gene_name_for_search])	#keyを利用してdictを回すための遺伝子名リスト
		if gene_name_for_search in line:
			annotation_dict_list.setdefault(gene_name_for_search, []).append(tuple(line))	#遺伝子名をkey,valueをannotationの行としてdictを作る(valueには複数の要素が入ってる)
	unique_gene_name_for_search_list = gene_id_list		#gene_idコードのリストを読み込みsearch対象とする
	unique_gene_name_for_search_list.sort()		#一応遺伝子番号順にソート
	print('input gene id number is', len(unique_gene_name_for_search_list))
	if desgin_pos == 'fs':
		feature_name_list = ['CDS','exon','miRNA']
	elif desgin_pos == 'ld':
		feature_name_list = ['gene']
	annotation_list = []
	for annotation_gene_name_line in unique_gene_name_for_search_list:	#keyを回すための遺伝子名
		flag = 'off'	#flagを初期化
		for feature_name_line in feature_name_list:	#featureリストをまわす
			if flag == 'on':	#優先順位の高いfeatureがリストから取れればあとはブレイク
				break
			try:
				for annotation_per_gene in list(annotation_dict_list[annotation_gene_name_line]):	#遺伝子名リストによりdictを回す
					[reference, source, feature, start, end, score, strand, frame, gene_name] = list(annotation_per_gene)
					if feature_name_line == feature:
						flag = 'on'	#featureと一致したらflagがonになる
						annotation_list.append([reference, source, feature, start, end, score, strand, frame, gene_name])
			except KeyError:
				continue
	print ('annotation list length is ' + str(len(annotation_list)))
	return annotation_list

def process_pam_func(inputted_pam):
	inputted_pam = inputted_pam
	pam_len = len(inputted_pam)
	fw_pam = inputted_pam.translate(str.maketrans('N', '.'))   #.GG
	rv_pam = fw_pam.translate(str.maketrans('ATGC', 'TACG'))[::-1]   #CC.
	return [fw_pam, rv_pam, pam_len]
