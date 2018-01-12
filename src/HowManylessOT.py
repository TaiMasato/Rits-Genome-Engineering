#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
if __name__ == '__main__':
	openfile = '../assets/selected_GC_gRNAs_Entire.txt'
	print openfile
	fewOT = 0
	manyOT = 0
	num =0
	for i in open(openfile, 'r'):
		num +=1
		i = i.rstrip()
		OTnum = re.split('[\t=]',i)[-1]
		if int(OTnum) == 0:
			fewOT += 1
		elif int(OTnum) != 0:
			manyOT += 1
			
	print 'OTが少ないgRNA'
	print fewOT
	print 'OTがあるgRNA'
	print manyOT
	print '比率'
	print (float(fewOT)/(float(fewOT)+float(manyOT)))*100,
	print '%'
	print num
		
		