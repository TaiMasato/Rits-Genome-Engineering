#-*- coding: utf-8 -*-
#!/usr/bin/env python
#日本語
import time
start = time.time()

def ChooseGCcontent(openfile):
	outside_range = 0
	inside_range = 0
	
	
	for gff_line in open(openfile , 'r'):
		gff_line = gff_line.rstrip()
		[reference , annotation , feature , start , end , score , strand , frame , attribute] = gff_line.split('\t')
		attributeitems = attribute.split(';')
		gRNAname = attributeitems[0]
		gRNAseq = attributeitems[1][:-3]  #PAM配列は除く
		Gcontent = gRNAseq.count('G')
		Ccontent = gRNAseq.count('C')
		GCcontent = (float(Gcontent + Ccontent)/20)*100
		
		if GCcontent < 25 or 75 < GCcontent:
			continue
		f.write(str(reference) + '\t' + str(annotation) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attribute) + '\n' )

		
if __name__ == '__main__':
	openfile = '../assets/UniquegRNAlist_withAnnotation.txt'
	outputfile = '../assets/selected_GC_gRNAs.txt'
	f = open(outputfile,'w')
	ChooseGCcontent(openfile)
	f.close
	print time.time() - start,
	print 'sec : Done'
	