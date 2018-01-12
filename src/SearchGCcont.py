#-*- coding: utf-8 -*-
#!/usr/bin/env python
#日本語

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
		#print str(gRNAname) + '\t' + str(gRNAseq) + '\t' +str(GCcontent) + '%'     #各gRNAのGC含量を見るときはここを解除(本数多いとどエラいことになるので注意)
		if int(GCcontent) < 25 or 75 < int(GCcontent):
			outside_range += 1
		elif 25 <= int(GCcontent) and int(GCcontent) <= 75:
			inside_range += 1
	
		"""
		if GCcontent < 25 and 75 < GCcontent:
			continue
		f.write(str(reference) + '\t' + str(annotation) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attribute) + '\n' )
		print gRNAname,
		print gRNAseq,
		print (float(Gcontent + Ccontent)/20)*100
		"""
	print '範囲外:',
	print outside_range
	print '範囲内:',
	print inside_range
	print 'GC含量が範囲内のgRNAが存在する割合',
	print float(inside_range/float(outside_range + inside_range))*100,
	print '%'
if __name__ == '__main__':
	openfile = '../assets/selected_gRNAs.txt'
	print openfile
	ChooseGCcontent(openfile)