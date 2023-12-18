import numpy as np
import pickle
import regex
import sys
import os


def gefa(new_fafile,old_fafile):
    with open(new_fafile,'w') as outputfile:
        outputfile.write('')
    with open(old_fafile,'r') as fafile:
        content = fafile.read()
    gene_list =content.split('>')
    del gene_list[0]
    for gene in gene_list:
        gene_name = gene.split('\n',1)[0]
        gene_sequence = gene.split('\n',1)[1].replace('\n','')
        gene_infomation = gene_name.split(' ')
        Transcrtpt = gene_infomation[0].split('.')[0]
        chr_name = 'chr'+gene_infomation[2].split(':')[2]
        start = gene_infomation[2].split(':')[3]
        end = gene_infomation[2].split(':')[4]
        strand = gene_infomation[2].split(':')[5]
        if strand == str(1):
            strand = '+'
        elif strand == str(-1):
            strand ='-'
        Gene = gene_infomation[3].split('.')[0].split(':')[1]
        Type = gene_infomation[5].split('.')[0].split(':')[1]
        if 'CHR' in chr_name or '.' in chr_name:
            continue
        if 'MT' in chr_name:
            continue
        if 'gene_symbol' not in gene_name:
            Symbol = 'no Symbol'
        else:
            Symbol = gene_infomation[6].split(':')[1]
        #打印转录本信息
        #print(Transcrtpt,Gene,Type,Symbol,chr_name,start,end,strand)
        with open(new_fafile,'a')as outputfile:
            outputfile.write('>' + Transcrtpt + ';' + Gene + ';' + Type + ';' + Symbol + ';' + chr_name + ';' + start + ';' + end + ';' + strand + '\n')
            outputfile.write(gene_sequence+'\n')

#读出染色体文件
def readchr(chrfile):
    with open(chrfile,'r')as inputfile:
        content = inputfile.read()
    chr_list = content.split('>')
    chr_dict = dict()
    del chr_list[0]
    for chr in chr_list:
        chr_name = chr.split('\n',1)[0]
        chr_sequence = chr.split('\n',1)[1].replace('\n','')
        if 'CHR' in chr_name or '.' in chr_name:
            continue
        if 'MT' in chr_name:
            continue
        chr_dict[chr_name] = chr_sequence
    return chr_dict

def readbed(bed_file):
    with open(bed_file,'r') as inputfile:
        content = inputfile.read()
    bed_list = content.split('>')
    bed_dict = {}
    del bed_list[0]
    for bed in bed_list:
        bed = bed.split('\n',1)[0]
        bed_information = bed.split(';')
        chr_name = bed_information[4]
        start = bed_information[5]
        end = bed_information[6]
        transcript = bed_information[0].split('.')[0]
        gene = bed_information[1]
        strand = bed_information[7]
        type = bed_information[2]
        symbol = bed_information[3]
        if strand == '+':
            promoter_start = int(start)-2000
            promoter_end = int(start)+2000
        if strand == '-':
            promoter_start = int(end)+2000
            promoter_end = int(end)-2000
        bed_dict[transcript] = {}
        bed_dict[transcript]['chr_name'] = chr_name
        bed_dict[transcript]['start'] = start
        bed_dict[transcript]['end'] = end
        bed_dict[transcript]['gene'] = gene
        bed_dict[transcript]['strand'] = strand
        bed_dict[transcript]['type'] = type
        bed_dict[transcript]['symbol'] = symbol
        bed_dict[transcript]['promoter_start'] = promoter_start
        bed_dict[transcript]['promoter_end'] = promoter_end
    return bed_dict

def read_fa(fa_file_name):
    fa_dict = {}
    fa_len_dict = {}
    with open(fa_file_name, 'r') as fa_file:
        contents = fa_file.read()
        fa_list = contents.split('>')
        for i in range(1, len(fa_list)):
            fa = fa_list[i]
            info, sequence = fa.split('\n', 1)
            info = info.split()[0] # deal with multi column information
            sequence = sequence.replace('\n', '')
            fa_dict[info] = sequence.upper()
            fa_len_dict[info] = len(sequence)
    return (fa_dict, fa_len_dict)

def reverse_complement(sequence):#返回反向互补序列
    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
    sequence = sequence.translate(trantab)
    return sequence[::-1]

def gepromoter(withPromoterfile,new_fafile,bed_dict,chr_dict):
    with open(withPromoterfile,'w') as outputfile:
        outputfile.write('')
    with open(new_fafile,'r')as inputfile:
        content = inputfile.read()
    gene_list = content.split('>')
    del gene_list[0]
    for gene in gene_list:
        bed = gene.split('\n', 1)[0]
        sequence = gene.split('\n',1)[1].replace('\n','')
        bed_information = bed.split(';')
        transcript = bed_information[0]
        chr_name = bed_dict[transcript]['chr_name']
        strand = bed_dict[transcript]['strand']
        start = bed_dict[transcript]['start']
        end = bed_dict[transcript]['end']
        promoter_start = bed_dict[transcript]['promoter_start']
        promoter_end = bed_dict[transcript]['promoter_end']
        if strand == '+' and chr_name in chr_dict:
            with open(withPromoterfile,'a') as outputfile:
                outputfile.write('>'+bed+'\n')
                outputfile.write(chr_dict[chr_name][promoter_start-1:int(start)].upper()+sequence+'\n')
        if strand == '-' and chr_name in chr_dict:
            with open(withPromoterfile,'a') as outputfile:
                outputfile.write('>'+bed+'\n')
                outputfile.write(reverse_complement(chr_dict[chr_name][int(end)-1:promoter_start]).upper()+sequence+'\n')
        if chr_name not in chr_dict:
            print(bed)

def geconverse(new_fafile,format):
    with open(new_fafile, 'r')as inputfile:
        content = inputfile.read()
    gene_list = content.split('>')
    gene_dict = {}
    del gene_list[0]
    for gene in gene_list:
        bed = gene.split('\n', 1)[0]
        bed_information = bed.split(';')
        transcript = bed_information[0]
        gene = bed_information[1]
        symbol = bed_information[3]
        if format == 0:
            if gene not in gene_dict:
                gene_dict[gene] = []
                gene_dict[gene].append(transcript)
            else:
                gene_dict[gene].append(transcript)
        if format == 1:
            if symbol not in gene_dict:
                gene_dict[symbol] = []
                gene_dict[symbol].append(transcript)
            else:
                gene_dict[symbol].append(transcript)
    return gene_dict


def parse_gtf(gtf_file_name, chrom_len_dict):
    is_gff = False
    if gtf_file_name[-4:].upper() == '.GFF' or gtf_file_name[-5:].upper() == '.GFF3':
        is_gff = True
        id_dict = {}
        with open(gtf_file_name, 'r') as gtf_file:
            for line in gtf_file:
                if line[0] == '#':
                    continue
                items = line.strip().split('\t')
                try:
                    id_ = regex.findall(ID_PATTERN, items[-1])[0]
                    parent = regex.findall(PARENT_ID_PATTERN, items[-1])[0]
                    id_dict[id_] = parent
                except:
                    pass

    exon_dict = {}  # {chrom: {gene: {start_end, ...}, ...}, ...}
    with open(gtf_file_name, 'r') as gtf_file:
        for line in gtf_file:
            if line[0] == '#':
                continue
            items = line.strip().split('\t')
            bio_type = items[2]
            if bio_type == 'exon':
                # print(items)
                if 'CHR' in items[0] or '.' in items[0]:
                    continue
                if 'MT' in items[0]:
                    continue
                if 'DNA' in items[0]:
                    continue
                chrom = 'chr' + items[0]
                start = items[3]
                end = items[4]
                if is_gff:
                    try:
                        id_ = regex.findall(ID_PATTERN, items[-1])[0]
                        gene = id_dict[id_dict[id_]]
                    except:  # for exon's parent is pseudogene, skip
                        continue
                else:
                    gene = regex.findall(GENE_ID_PATTERN, items[-1])[0]
                if chrom not in exon_dict:
                    exon_dict[chrom] = {}
                if gene not in exon_dict[chrom]:
                    exon_dict[chrom][gene] = set()
                exon_dict[chrom][gene].add(start + '\t' + end)

    location_dict = {}  # {chrom: [[start, end, bio_type, gene], ...], ...}
    gene_dict = {}  # {chrom: [[start, end, gene], ...], ...}
    for chrom in exon_dict:
        chrom = chrom
        location_dict[chrom] = []
        gene_dict[chrom] = []

        # handle exon and intron
        for gene in exon_dict[chrom]:
            exon_set = exon_dict[chrom][gene]
            exon_list = [list(map(int, exon.split('\t'))) for exon in exon_set]

            # get gene block
            gene_start = min(set(sum(exon_list, [])))
            gene_end = max(set(sum(exon_list, [])))
            gene_dict[chrom].append([gene_start, gene_end, gene])

            # remove overlap
            exon_list_remove_overlap = []
            exon_list.sort()
            for i in range(len(exon_list) - 1):
                if exon_list[i + 1][0] < exon_list[i][1]:
                    exon_list[i + 1][0] = exon_list[i][0]
                    exon_list[i + 1][1] = max(exon_list[i][1], exon_list[i + 1][1])
                else:
                    exon_list_remove_overlap.append(exon_list[i])
            exon_list_remove_overlap.append(exon_list[-1])

            # add exon
            for exon in exon_list_remove_overlap:
                exon.extend(['exon', gene])
                location_dict[chrom].append(exon)
            # get intron
            for i in range(len(exon_list_remove_overlap) - 1):
                intron = [exon_list_remove_overlap[i][1] + 1, exon_list_remove_overlap[i + 1][0] - 1, 'intron', gene]

                location_dict[chrom].append(intron)

        # handle intergenic
        gene_list = gene_dict[chrom]
        # add head and tail
        gene_list.append([0, 0, 'chrom_head'])
        chrom_len = chrom_len_dict[chrom]
        gene_list.append([chrom_len + 1, chrom_len + 1, 'chrom_tail'])
        # remove overlap
        gene_list.sort()

        gene_list_remove_overlap = []
        for i in range(len(gene_list) - 1):
            if gene_list[i + 1][0] < gene_list[i][1]:
                gene_list[i + 1][0] = gene_list[i][0]
                gene_list[i + 1][1] = max(gene_list[i][1], gene_list[i + 1][1])
                gene_list[i + 1][2] = gene_list[i][2] + '_' + gene_list[i + 1][2]
            else:
                gene_list_remove_overlap.append(gene_list[i])
        gene_list_remove_overlap.append(gene_list[-1])

        # get intergenic
        for i in range(len(gene_list_remove_overlap) - 1):
            gene = gene_list_remove_overlap[i][2] + ' & ' + gene_list_remove_overlap[i + 1][2]
            intergenic = [gene_list_remove_overlap[i][1] + 1, gene_list_remove_overlap[i + 1][0] - 1, 'intergenic',
                          gene]

            location_dict[chrom].append(intergenic)

    for chrom in location_dict:
        location_dict[chrom].sort()

    # print(location_dict['Chr1'])
    return location_dict


def binary_search(query, subject): # query [start, end], subject [[start, end, bio_type, gene], ...]
    query_start, query_end = query
    low = 0
    high = len(subject) - 1
    mid = 0
    while high >= low:
        mid = (low + high) // 2
        #print(mid)
        candidate = subject[mid]
        #print(candidate)
        candidate_start, candidate_end, bio_type = candidate[:3]
        if (candidate_start <= query_start and query_end <= candidate_end): # query total in candidate [candidate_start, query_start, query_end, candidate_end]
            # print(1)
            return candidate # [start, end, bio_type, gene]
        elif (query_start < candidate_start and query_end >= candidate_start and query_end <= candidate_end): # query_ritht overlap candidate_left
            if bio_type == 'exon': # if overlap, return exon
                # print(2)
                return candidate
            else:
                # print(3)
                return  subject[mid - 1]
        elif (candidate_start <= query_start and query_start <= candidate_end and candidate_end < query_end): # query_left overlap candidate_right
            if bio_type == 'exon': # if overlap, return exon
                # print(4)
                return candidate
            else:
                # print(5)
                return  subject[mid + 1]
        elif (query_start < candidate_start and candidate_end < query_end): # query_left overlap candidate_right
            if bio_type == 'exon': # if overlap, return exon
                # print(6)
                return candidate
            else:
                # print(7)
                return  subject[mid + 1]
        elif (query_end < candidate_start): # continue searching on the left
            high = mid - 1
            # print('left')
        elif (query_start > candidate_end): # continue searching on the right
            low = mid + 1
            # print('right')
        else:
            print('binary search error')
            # print(candidate)
    range_start = mid - 1000
    range_end = mid + 1000
    if range_start < 0:
        range_start = 0
    if range_end >= len(subject):
         range_end = len(subject) - 1
    for mid in range(range_start, range_end):
        candidate = subject[mid]
        # print(candidate)
        candidate_start, candidate_end, bio_type = candidate[:3]
        if (candidate_start <= query_start and query_end <= candidate_end): # query total in candidate [candidate_start, query_start, query_end, candidate_end]
            # print(1)
            return candidate # [start, end, bio_type, gene]
        elif (query_start < candidate_start and query_end >= candidate_start and query_end <= candidate_end): # query_ritht overlap candidate_left
            if bio_type == 'exon': # if overlap, return exon
                # print(2)
                return candidate
            else:
                # print(3)
                return  subject[mid - 1]
        elif (candidate_start <= query_start and query_start <= candidate_end and candidate_end < query_end): # query_left overlap candidate_right
            if bio_type == 'exon': # if overlap, return exon
                # print(4)
                return candidate
            else:
                # print(5)
                return  subject[mid + 1]
        elif (query_start < candidate_start and candidate_end < query_end): # query_left overlap candidate_right
            if bio_type == 'exon': # if overlap, return exon
                # print(6)
                return candidate
            else:
                # print(7)
                return  subject[mid + 1]
    return '' # failure

GENE_ID_PATTERN = 'gene_id "([^;]+)";' # for GTF
ID_PATTERN = 'ID=([^;]+)' # for GFF
PARENT_ID_PATTERN = 'Parent=([^;]+)' # for GFF

#这里修改文件
# OLD_FA_FILE = 'Annotation/Danio_rerio.GRCz10.ncrna.fa' # 需要的原始文件1，非编码基因序列。
# NEW_FA_FILE = 'Annotation/Zebrafish_ncrna.fa'  # 生成的中间文件1，非编码基因序列整合版。
# BED_FA = 'Annotation/Zebrafish_ncrna.fa'  # 生成的中间文件1，非编码基因序列整合版。
# CHR_FA = 'fasta/Zebrafish.fa'  # 需要的原始文件2，基因组序列。
# WITH_PROMOTER_FILE = 'Annotation/Zebrafish_ncRNA_withPromoter.fa' # 生成的中间文件2，带有启动子序列的序列文件。
# GTF_FILE_NAME = 'Annotation/Danio_rerio.GRCz10.80.gtf'  # 需要的原始文件3，ensemble下载的gtf注释文件。
# BED_DICT_FILE = 'dict/Zebrafish_bed_dict.pkl'  # 生成的中间文件3，基因的bed注释。
# CHR_DICT_FILE = 'dict/Zebrafish_chr_dict.pkl'  # 生成的中间文件4，染色体的dict文件。
# LOCATION_DICT_FILE = 'dict/Zebrafish_location_dict.pkl' #生成的中间文件5，位置区域的dict文件。

OLD_FA_FILE = sys.argv[1] # 需要的原始文件1，非编码基因序列。
GTF_FILE_NAME = sys.argv[2] # 需要的原始文件2，ensemble下载的gtf注释文件
CHR_FA = sys.argv[3] # 需要的原始文件3，基因组序列


#输出文件，统一命名
NEW_FA_FILE = 'Annotation/ncrna.fa'  # 输出文件1，非编码基因序列整合版。
BED_FA = 'Annotation/ncrna.fa'  # 输出文件1，非编码基因序列整合版。
WITH_PROMOTER_FILE = 'Annotation/ncRNA_withPromoter.fa' # 输出文件2，带有启动子序列的序列文件。
BED_DICT_FILE = 'dict/bed_dict.pkl'  # 输出文件3，基因的bed注释。
CHR_DICT_FILE = 'dict/chr_dict.pkl'  # 输出文件4，染色体的dict文件。
LOCATION_DICT_FILE = 'dict/location_dict.pkl' #输出文件5，位置区域的dict文件。

if __name__ == '__main__':

    #第一步生成新的fa文件,只用运行一次
    gefa(NEW_FA_FILE,OLD_FA_FILE)

    #第二步保存中间文件，只用运行一次
    chr_dict = readchr(CHR_FA)
    f_save = open(CHR_DICT_FILE, 'wb')
    pickle.dump(chr_dict, f_save)
    f_save.close()
    bed_dict = readbed(BED_FA)
    f_save = open(BED_DICT_FILE, 'wb')
    pickle.dump(bed_dict, f_save)
    f_save.close()
    fa_dict , fa_len_dict = read_fa(CHR_FA)
    location_dict = parse_gtf(GTF_FILE_NAME, fa_len_dict)
    f_save = open(LOCATION_DICT_FILE, 'wb')
    pickle.dump(location_dict, f_save)
    f_save.close()
    #第三步读取中间的chr_dict和bed_dict,location_dict
    #加载文件
    chr_read = open(CHR_DICT_FILE, 'rb')
    chr_dict = pickle.load(chr_read)
    bed_read = open(BED_DICT_FILE, 'rb')
    bed_dict = pickle.load(bed_read)
    f_read = open(LOCATION_DICT_FILE, 'rb')
    location_dict = pickle.load(f_read)


    #第四步 生成promoter文件
    gepromoter(WITH_PROMOTER_FILE,NEW_FA_FILE,bed_dict,chr_dict)

