#-*- coding: utf-8 -*-
#!/usr/bin/env python
#日本語

#殴り書きです

import time
import re
start = time.time()


def eachOTcount(pickupfile):
	OTlist = []
	pregene = ""
	offtarget = ""
	for fileline in open(pickupfile , 'r'):
		fileline =  fileline.rstrip()
		[reference , source , feature , start , end , score , strand , frame , attribute] =  fileline.split('\t')
		attribute = re.split('[;=]' , attribute)
		if pregene == attribute[0]:
			continue
		elif pregene != attribute[0] and attribute[0][-1] == "u":
			pregene = attribute[0]
			offtarget = str(attribute[12])
		elif pregene != attribute[0] and attribute[0][-1] == "d":
			pregene = attribute[0]
			OTlist.append([offtarget+str(attribute[12])])
	return OTlist

def aaa(OTlist):
	withOT = 0
	withoutOT = 0

	for OTnum in OTlist:
		if int(OTnum[0]) != 00:
			withOT += 1
		elif int(OTnum[0]) == 00:
			withoutOT += 1
	print 'OTなしのgRNAの確率',
	print (float(withoutOT)*100 / (float(withOT) + float(withoutOT))),
	print '%'


if __name__ == '__main__':
	pickupfile = '../assets/selected_gRNAs_Mp.txt'
	OTlist = eachOTcount(pickupfile)
	aaa(OTlist)
