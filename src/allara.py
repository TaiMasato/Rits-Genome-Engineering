#!/usr/bin/env python
# -*- coding: utf-8 -*-
#スプライシングバリアント考慮有り
import time
import itertools
import string
import re
import commands

main_genome_fasta = '../assets/sca119.txt'
mapoly_gff3 = '../assets/Mpolymorphav3.1.gene.gff3'
#mapoly_gff3 = '../assets/made_gff.txt'
output_file = '../assets/allgRNAList_GFFformat.txt'
output_file_for_bowtie = '../assets/allgRNAList_for_bowtie.txt'
output_file_fastaname = '../assets/fastaname.txt'
output_file_fastaseq = '../assets/fastaseq.txt'
get_PAM = 'NGG'
OneBeforeGeneName = 'genename'
allgrna_in_gene_no=1

def MainGenomeFileOpen(GenomeFastaFile):
	'''ゼニゴケのmainゲノムfastaファイルから遺伝子名、seqだけのlistを作る'''
	GenomeFastaNameList = []
	GenomeFastaFragSeq = []
	GenomeFastaEachSeq= []
	GenomeFastaSeqList = []
	
	Genome_fasta_file = open(GenomeFastaFile)   #(1)
	Genomefasta__line = Genome_fasta_file.readline().rstrip()   #(2)
	if Genomefasta__line.startswith(">"):   #(3)
		GenomeFastaNameList.append(Genomefasta__line)
	for Genomefasta__line in Genome_fasta_file:
	 	Genomefasta__line = Genomefasta__line.strip()   #(4)
	 	if '>' in Genomefasta__line:   #(5)
	 		GenomeFastaEachSeq.append(''.join(GenomeFastaFragSeq))   #(6)
	 		GenomeFastaSeqList.extend(GenomeFastaEachSeq)
	 		GenomeFastaEachSeq = []
	 		GenomeFastaFragSeq = []	#(7)
	 		GenomeFastaNameList.append(Genomefasta__line)
	 	else:   #(8)
	 		GenomeFastaFragSeq.append(Genomefasta__line)   #(9)
	GenomeFastaEachSeq.append(''.join(GenomeFastaFragSeq))   #(10)
	GenomeFastaSeqList.extend(GenomeFastaEachSeq)
	GenomeFastaEachSeq = []
	
	for i in xrange(len(GenomeFastaNameList)):
		f3.write(str(GenomeFastaNameList[i]) + '\n')
		f4.write(str(GenomeFastaSeqList[i]) + '\n')
	
	
	return [GenomeFastaNameList , GenomeFastaSeqList]
		
def Mapolygff3Open(gff3File):
	gff3List = []    #gff3をそのままリスト化したもの
	referenceList = []   #gff3からreference要素のみを集めてリスト化したもの
	preReferenceName = 'reference'
	referenceCounter = 0
	
	for gff_line in open(gff3File, 'r'): 
		if "##" in gff_line: 
			continue
		gff_line = gff_line.rstrip()
		[reference, source, feature, start, end, score, strand, frame, attr] = gff_line.split('\t')
		
		if feature == "CDS":
			if reference != preReferenceName:
				referenceList.append(reference)
				referenceCounter += 1
				preReferenceName = reference
			attr = attr.split('.')
			attr = attr[0][3:] + '.' + attr[1]
			gff3List.append([reference, source, feature, start, end, score, strand, frame, attr])
	return [gff3List , referenceCounter , referenceList]	
	
def ProcessPAM(get_pam):
	pam_len = len(get_pam)
	fw_pam=get_pam.translate(string.maketrans('N', '.'))   #.GG
	rv_pam = fw_pam.translate(string.maketrans('ATGC', 'TACG'))[::-1]   #CC.
	return [fw_pam , rv_pam , pam_len]

def SearchgRNA(gff_reference , anotation , feature , cds_start , cds_end , score , strand , frame , gene_name , CDSeq , PAM , PamLen ):
	global OneBeforeGeneName
	global allgrna_in_gene_no
	
	match = re.compile(PAM)
	i = 0
	while i >= 0:
		matching = match.search(CDSeq, i)
		if not matching:
			break
		#print matching.start()
		if OneBeforeGeneName != gene_name:
			allgrna_in_gene_no = 1
			OneBeforeGeneName = gene_name
		PamPosition = matching.start()
		FwSpacerStart = PamPosition - 20
		RvSpacerStart = PamPosition + 3
		if PAM == '.GG':
			grna_start_position = int(cds_start) + int (FwSpacerStart)
			gRNA = CDSeq[FwSpacerStart : matching.start()]
			if 'TTTT' not in gRNA and len(gRNA) >= 20:
				#print gRNA
				f1.write(str(gff_reference) + '\t' + str(anotation) + '\t' + str(feature) + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+20) + '\t' + '.' + '\t' + str(strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(strand) + '=' + str(grna_start_position) + '=' + str(allgrna_in_gene_no) + ';' + str(gRNA) + '\n' )
				f2.write('>' + str(gene_name) + '=' + str(strand) + '=' + str(grna_start_position) + '=' + str(allgrna_in_gene_no) + '\n' + str(gRNA) + '\n' )
			allgrna_in_gene_no += 1
		
		if PAM == 'CC.':
			grna_start_position = int(cds_start) + int (RvSpacerStart)
			gRNA = CDSeq[RvSpacerStart : RvSpacerStart + 20]
			if 'AAAA' not in gRNA and len(gRNA) >= 20:
				gRNA = gRNA.translate(string.maketrans('ATGC', 'TACG'))[::-1]
				#print gRNA
				f1.write(str(gff_reference) + '\t' + str(anotation) + '\t' + str(feature) + '\t' + str(grna_start_position) + '\t' + str(int(grna_start_position)+20) + '\t' + '.' + '\t' + str(strand) + '\t' + '.' + '\t' + str(gene_name) + '=' + str(strand) + '=' + str(grna_start_position) + '=' + str(allgrna_in_gene_no) + ';' + str(gRNA) + '\n' )
				f2.write('>' + str(gene_name) + '=' + str(strand) + '=' + str(grna_start_position) + '=' + str(allgrna_in_gene_no) + '\n' + str(gRNA) + '\n' )
			allgrna_in_gene_no += 1
		i = matching.start() + 1

def MakegRNA(ReferenceNameList , GenomeSeqList , gff3List , FwPam , RvPam , Pam_Len ):
	for x , y in itertools.product(xrange(len(gff3List)), xrange(len(GenomeSeqList))):
		[gff_reference , anotation , feature , start , end , score , strand , frame , gene_name] = [gff3List[x][0] , gff3List[x][1] , gff3List[x][2] , gff3List[x][3] , gff3List[x][4] , gff3List[x][5] , gff3List[x][6] , gff3List[x][7] , gff3List[x][8]]
		main_genome_reference = ReferenceNameList[y][1:]
		main_genome_seq = GenomeSeqList[y]
		if gff_reference != main_genome_reference:
			continue
		cdseq = main_genome_seq[int(start)-1:int(end)]
		#print cdseq
		SearchgRNA(gff_reference , anotation , feature , start , end , score , strand , frame , gene_name , cdseq , FwPam , Pam_Len )
		SearchgRNA(gff_reference , anotation , feature , start , end , score , strand , frame , gene_name , cdseq , RvPam ,  Pam_Len )

		
def BowtieCheck():
   	f = open("../assets/BowtieResultFile.txt","w")
   	bowtie_result = []
   	bowtie_result = (commands.getoutput('bowtie -a -v 3 ../../../../..//usr/local/share/bowtie/indexes/Mapoly --suppress 6,7 -f ../assets/allgRNAList_for_bowtie.txt'))
   	#bowtie_result = (commands.getoutput('bowtie -a -v 3 ../bowtie-1.2.1.1/indexes/mapoly --suppress 6,7 -f ../assets/allgRNAList_for_bowtie.txt'))	#Linux用
   	bowtie_result_list = bowtie_result
   	f.write(bowtie_result)
   	f.close()
   	
def BowtieResultFixer():
	f = open("../assets/BowtieResultFixedFile.txt","w")
	for line in open("../assets/BowtieResultFile.txt","r"):
		if line.startswith('#'):
			break
		f.write(line)
	f.close()
	
if __name__ == '__main__':
	'''各種関数'''
	TimeMeasurement = time.time()

	f1 = open(output_file, 'w')
	f2 = open(output_file_for_bowtie, 'w')
	f3 = open(output_file_fastaname , 'w')
	f4 = open(output_file_fastaseq , 'w')
	
	GenomeFastaNameList , GenomeFastaSeqList = MainGenomeFileOpen(main_genome_fasta)
	gff3List , referenceCounter , referenceList = Mapolygff3Open(mapoly_gff3)
	fw_pam , rv_pam , pam_len = ProcessPAM(get_PAM)
	
	f1.write ('# reference'+ '\t'  + 'source'+ '\t' + 'feature' + '\t'  + 'start' + '\t'+ 'end' + '\t' + 'score' + '\t' + 'strand' + '\t' + 'frame' +'\t' + '色々' + '\n')
	MakegRNA(GenomeFastaNameList , GenomeFastaSeqList , gff3List , fw_pam , rv_pam , pam_len)
	f1.close
	f2.close
	f3.close
	f4.close
	print 'Making gRNA Finished'
	print time.time() - TimeMeasurement,
	print 'sec'
	#BowtieCheck()
	#BowtieResultFixer()
	print 'Bowtie check finished'
	print time.time() - TimeMeasurement,
	print 'sec'
	
	
	
