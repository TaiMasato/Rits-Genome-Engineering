
fastaname = '../assets/fastaname.txt'
fastaseq = '../assets/fastaseq.txt'

name = [i.rstrip().split('\t') for i in open(fastaname,'r')]
seq = [j.rstrip().split('\t') for j in open(fastaseq, 'r')]





for i in xrange(len(name)):
	if '>scaffold_119' == name[i][0]:
		print seq[i]