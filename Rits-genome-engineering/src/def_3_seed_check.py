import re
import time

def seed_check_func(output_pam_check, output_file):
	'''
	bowtieによって出力されたファイルをリスト化する。
	行と列を入れ替えたものを別のリストに格納し、各遺伝子に幾つのofftargetがあるかをカウント。
	offtargetの数を元のリストに新要素として格納、拡張する。
	'''
	TimeMeasurement = time.time()       #処理時間の計算開始
	BowtieResult_with_MMposition_line = []
	BowtieResult_ItemList = []
	BowtieResult_with_MMposition_List = []
	tr_BowtieResult_with_MMposition_List = []
	BowtieResult_with_MMcount_List = []
	seedarea = 12
	seed_num = 0
	nonseed_num = 0
	f = open(output_file,"w")
	for BowtieResultFile_line in open(output_pam_check, 'r'):   #ファイルの展開
		BowtieResultFile_line = BowtieResultFile_line.strip()
		BowtieResultFile_line = BowtieResultFile_line.split('\t')
		#[GeneName, strand, scaffold, position, sequence,mismatch]
		try:
			[GeneName, strand, scaffold, position, sequence,mismatch] = BowtieResultFile_line
		except ValueError:
			[GeneName, strand, scaffold, position, sequence] = BowtieResultFile_line
			mismatch =''
		sequence = sequence.strip()
		mismatch = mismatch.strip()
		MutationCounter = mismatch.count(':')
		'''
		Perfect matchのケース
		'''
		if MutationCounter == 0:
			BowtieResult_with_MMposition_line.extend([GeneName, strand, scaffold, position, sequence, mismatch,str(MutationCounter).zfill(2),str(seed_num).zfill(2),str(nonseed_num).zfill(2)])   #整頓したものを再リスト化(１次元)
			for element in BowtieResult_with_MMposition_line:
				f.write(str(element)+'\t')
			f.write(str('\n'))
			#BowtieResult_with_MMposition_List.append(BowtieResult_with_MMposition_line)
			 #再リスト化(２次元)
			BowtieResult_with_MMposition_line=[]
			continue
		'''
		Mutationが一つ以上のケース
		'''
		splitted_mismatch = re.split('[,:]',mismatch)
		for x in range(MutationCounter):
			if int(splitted_mismatch[2*x]) < seedarea:
				seed_num += 1
				continue
			nonseed_num += 1
		BowtieResult_with_MMposition_line.extend([GeneName, strand, scaffold, position, sequence,mismatch,str(MutationCounter).zfill(2),str(seed_num).zfill(2),str(nonseed_num).zfill(2)])   #整頓したものを再リスト化(１次元)
		for element in BowtieResult_with_MMposition_line:
			f.write(str(element)+'\t')
		f.write(str('\n'))
		#BowtieResult_with_MMposition_List.append(BowtieResult_with_MMposition_line)
		BowtieResult_with_MMposition_line=[]
		seed_num = 0
		nonseed_num = 0
	#csvWriter = csv.writer(f)
	print ('Miss Match Check Done!')
	print ((time.time() - TimeMeasurement), 'sec')
