import re
import time
def gene_id_count(desgin_pos, gene_id_list, species_num):
	gene_name_list = []
	for gene_id_line in gene_id_list:
		if desgin_pos == 'fs':
			if species_num == '1':
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}\.\d', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}\.\d', gene_id_line)	#ゼニゴケY染色体の正
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d\.\d)', gene_id_line)
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '4':		#ポプラのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
		elif desgin_pos == 'ld':
			if species_num == '1':
				gene_id = re.search('Mapoly[\d]{4}s[\d]{4}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
				if gene_id == None:
					gene_id = re.search('Mapoly\w_\w[\d]{4}', gene_id_line)	#ゼニゴケY染色体の正
			elif species_num == '0':		#ナズナのgene_idcodeを取ってくる　総数:59516
				gene_id = re.search('(\w\w\w\w\d\d\d\d\d)', gene_id_line )
			elif species_num == '2':		#タバコのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Niben101Scf[\d]{5}g[\d]{5}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '3':		#トマトのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Solyc[\d]{2}g[\d]{6}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
			elif species_num == '4':		#ポプラのgene_idcodeを取ってくる	総数:24489
				gene_id = re.search('Potri.[\d]{3}G[\d]{6}', gene_id_line)	#ゼニゴケ常染色体とX染色体の正規表現
		if gene_id:
			gene_name = gene_id.group()
			if desgin_pos == 'fs':
				gene_name_list.extend([gene_name])
			elif desgin_pos == 'ld':
				gene_name_list.extend([gene_name + '=u'])
				gene_name_list.extend([gene_name + '=d'])
	gene_name_list = list(set(gene_name_list))	#重複削除
	gene_name_list.sort()		#遺伝子名順にソート
	return gene_name_list

def make_gRNA_dict_list(desgin_pos, gRNA_file):
	gRNA_dict_list = {}
	grna_list = [j.rstrip().split() for j in open(gRNA_file, 'r')]
	grna_list.sort(key=lambda x:x[8].split('=')[2])		#cut_positionでソート
	grna_list.sort(key=lambda x:x[8].split('=')[0])		#遺伝子名でソート
	for grna_line in grna_list:
		if desgin_pos == 'fs':
			gRNA_dict_list.setdefault(grna_line[8].split('=')[0], []).append(tuple(grna_line))
		elif desgin_pos == 'ld':
			gRNA_dict_list.setdefault(grna_line[8].split('=')[0] + '=' + grna_line[8].split('=')[1], []).append(tuple(grna_line))
	return gRNA_dict_list

def select_gRNAs_to_final(gRNA_dict_list, gene_id_list, gRNA_taking_num, f):
	cannot_obtain_the_required_number_of_gRNAs = []
	key_error_list = []
	gRNA_taking_num = int(gRNA_taking_num)
	for gene_id_code in gene_id_list:
		gRNA_list_per_gene = []
		try:
			for gRNA_dict_line in list(gRNA_dict_list[gene_id_code]):
				gRNA_list_per_gene.extend([list(gRNA_dict_line)])
				gRNA_list_per_gene.sort(key=lambda x:x[8].split('=')[-1])
			gRNA_counter = len(gRNA_list_per_gene)
			if gRNA_taking_num <= gRNA_counter:
				for x in range(int(gRNA_taking_num)) :
					[reference, source, feature, start, end, score, strand, frame, attr] = gRNA_list_per_gene[x]
					f.write(str(reference) + '\t' + str(source) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attr) + '\n')
			else:
				for x in range(int(gRNA_counter)) :
					[reference, source, feature, start, end, score, strand, frame, attr] = gRNA_list_per_gene[x]
					f.write(str(reference) + '\t' + str(source) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attr) + '\n')


			"""
			try:
				for x in range(int(gRNA_taking_num)) :
					[reference, source, feature, start, end, score, strand, frame, attr] = gRNA_list_per_gene[x]
					f.write(str(reference) + '\t' + str(source) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attr) + '\n')
			except IndexError:
				cannot_obtain_the_required_number_of_gRNAs.extend([gene_id_code])
				for x in range(len(gRNA_list_per_gene)) :
					[reference, source, feature, start, end, score, strand, frame, attr] = gRNA_list_per_gene[x]
					f.write(str(reference) + '\t' + str(source) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attr) + '\n')
			"""
		except KeyError:
			key_error_list.extend([gene_id_code])
	print ('# Can not obtain the required number of gRNAs list')
	print (cannot_obtain_the_required_number_of_gRNAs)
	print ('# Key error gene list')
	print (key_error_list)

def pick_out_gRNA_func(desgin_pos, gene_id_list, species_num, openfile, outputfile, extracted_num):
	TimeMeasurement = time.time()       #処理時間の計算開始
	gene_name_list = gene_id_count(desgin_pos, gene_id_list, species_num)
	gRNA_dict_list = make_gRNA_dict_list(desgin_pos, openfile)
	f = open(outputfile , 'w')
	select_gRNAs_to_final(gRNA_dict_list, gene_name_list, extracted_num, f)
	f.close()
	print ('Extract <n>gRNAs/Gene Done!')
	print ((time.time() - TimeMeasurement), 'sec')
