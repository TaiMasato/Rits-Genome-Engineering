import time


def open_gRNA_file(output_unique_gRNA_list_with_annotation):
    gRNA_list = [gff_line.rstrip().split('\t') for gff_line in open(output_unique_gRNA_list_with_annotation , 'r')]
    print (len(gRNA_list), 'gRNA number before processing')
    return gRNA_list
def restrict_GC_contents(gRNA_list):
    restrict_GC_gRNA_list = []
    for gff_line in gRNA_list:
        [reference, annotation, feature, start, end, score, strand, frame, attribute] = gff_line
        attributeitems = attribute.split(';')
        gRNAname = attributeitems[0]
        gRNAseq = attributeitems[1]
        Gcontent = gRNAseq.count('G')
        Ccontent = gRNAseq.count('C')
        GCcontent = (float(Gcontent + Ccontent)/len(gRNAseq))*100
        if GCcontent < 25 or 75 < GCcontent:		#GC含量が25~75%のgRNAを選出
            continue
        restrict_GC_gRNA_list.append(gff_line)
    return restrict_GC_gRNA_list

def restrict_poly_T(gRNA_list):
    restrict_poly_T_gRNA_list = []
    for gff_line in gRNA_list:
        [reference, annotation, feature, start, end, score, strand, frame, attribute] = gff_line
        attributeitems = attribute.split(';')
        gRNAname = attributeitems[0]
        gRNAseq = attributeitems[1]
        if 'TTTT' in gRNAseq:		#poly_Tを含むgRNAを削除
            continue
        restrict_poly_T_gRNA_list.append(gff_line)
    return restrict_poly_T_gRNA_list

def restriction_enzyme_site():
    pass

def output_restrict_gRNA(gRNA_list, outputfile):
    f = open(outputfile,'w')
    for gff_line in gRNA_list:
        [reference, annotation, feature, start, end, score, strand, frame, attribute] = gff_line
        f.write(str(reference) + '\t' + str(annotation) + '\t' + str(feature) + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + str(strand) + '\t' + str(frame) + '\t' + str(attribute) + '\n' )
    f.close()

def restrict_restrict_gRNA_func(restrict_gc, restrict_poly_t, output_unique_gRNA_list_with_annotation, outputfile):
    TimeMeasurement = time.time()       #処理時間の計算開始
    gRNA_list = open_gRNA_file(output_unique_gRNA_list_with_annotation)
    before_gRNA_len = len(gRNA_list)
    if restrict_gc == 'y':
        gRNA_list = restrict_GC_contents(gRNA_list)
        print (before_gRNA_len - len(gRNA_list), 'gRNAs was removed (GC content)')
        before_gRNA_len = len(gRNA_list)
    if restrict_poly_t == 'y':
        gRNA_list = restrict_poly_T(gRNA_list)
        print (before_gRNA_len - len(gRNA_list), 'gRNAs was removed (poly_T)')
#    if restrict_enzyme == 'y':
#        restriction_enzyme_site()
    output_restrict_gRNA(gRNA_list, outputfile)

    print ('Restricted gRNA List Done!')
    print ((time.time() - TimeMeasurement), 'sec')
