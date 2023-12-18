# coding:utf-8
import os, sys
import pickle
import regex
import numpy as np
# import RNA
import json
import re

import requests
from bs4 import BeautifulSoup
from requests.cookies import RequestsCookieJar

GENE_ID_PATTERN = 'gene_id "([^;]+)";'  # for GTF
ID_PATTERN = 'ID=([^;]+)'  # for GFF
PARENT_ID_PATTERN = 'Parent=([^;]+)'  # for GFF


# deal with promoterfa file
def read_promoterfa(filename, distance):
    promoterfa_dict = {}
    with open(filename, 'r') as fa_file:
        content = fa_file.read()
    fa_list = content.split('>')
    for i in range(1, len(fa_list)):
        # print(fa_list[i])
        info = fa_list[i]
        info_list = info.split('\n')[0].split(';')
        transcript = info_list[0]
        gene = info_list[1]
        type = info_list[2]
        symbol = info_list[3]
        chr_name = info_list[4]
        start = info_list[5]
        end = info_list[6]
        strand = info_list[7]
        sequence = info.split('\n', 1)[1].replace('\n', '')
        promoterfa_dict[transcript] = {}
        promoterfa_dict[transcript]['gene'] = gene
        promoterfa_dict[transcript]['type'] = type
        promoterfa_dict[transcript]['symbol'] = symbol
        promoterfa_dict[transcript]['chr_name'] = chr_name
        promoterfa_dict[transcript]['start'] = start
        promoterfa_dict[transcript]['end'] = end
        promoterfa_dict[transcript]['strand'] = strand
        promoterfa_dict[transcript]['sequence'] = sequence[2000 - distance:]

    return promoterfa_dict


# deal with fa file
def read_fa(filename):
    fa_dict = {}
    with open(filename, 'r') as fa_file:
        content = fa_file.read()
    fa_list = content.split('>')
    for i in range(1, len(fa_list)):
        # print(fa_list[i])
        info = fa_list[i]
        info_list = info.split('\n')[0].split(';')
        transcript = info_list[0]
        gene = info_list[1]
        type = info_list[2]
        symbol = info_list[3]
        chr_name = info_list[4]
        start = info_list[5]
        end = info_list[6]
        strand = info_list[7]
        sequence = info.split('\n', 1)[1].replace('\n', '')
        fa_dict[transcript] = {}
        fa_dict[transcript]['gene'] = gene
        fa_dict[transcript]['biotype'] = type
        fa_dict[transcript]['gene_symbol'] = symbol
        fa_dict[transcript]['chr_name'] = chr_name
        fa_dict[transcript]['start'] = start
        fa_dict[transcript]['end'] = end
        fa_dict[transcript]['strand'] = strand
        fa_dict[transcript]['sequence'] = sequence
    return fa_dict


def geconverse(new_fafile, format):
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


def save_obj(obj, name):
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# deal with sg_txt
def read_sg(file_name):
    length = 21
    with open(file_name, 'r') as sg_file:
        content = sg_file.read()
    sg_dict = {}
    gene_list = content.split('>')
    del gene_list[0]
    count = 1
    for gene in gene_list:
        gene_name = gene.split('\n', 1)[0].split(';')[0]
        gene_sequence = gene.split('\n', 1)[1].replace('\n', '')
        if gene_name not in sg_dict:
            sequence_num = 'sequence' + str(count)
            sg_dict[gene_name] = {}
            sg_dict[gene_name]['sequence'] = gene_sequence
            sg_dict[gene_name]['sequence_num'] = sequence_num
            count += 1

    return sg_dict


def Designsg(sequence):
    sg_dict = {}
    length = 21
    gene_sequence = sequence
    gene_sequence2 = reverse_complement(gene_sequence)
    match_objs = regex.finditer('[ATCG]{{{0}}}[AG]G'.format(length), gene_sequence)
    match2_objs = regex.finditer('[ATCG]{{{0}}}[AG]G'.format(length), gene_sequence2)
    for match_obj in match_objs:
        match_start, match_end = match_obj.span()
        sg_info = (match_obj.group(), (match_start, match_end, '+'))
        if match_obj.group() not in sg_dict:
            sg_dict[match_obj.group()] = [sg_info]
        else:
            sg_dict[match_obj.group()].append(sg_info)
    for match_obj2 in match2_objs:
        match_start, match_end = match_obj2.span()
        sg_info = (match_obj2.group(), (len(gene_sequence) - match_end, len(gene_sequence) - match_start, '-'))
        if match_obj2.group() not in sg_dict:
            sg_dict[match_obj2.group()] = [sg_info]
        else:
            sg_dict[match_obj2.group()].append(sg_info)
    # print(sg_dict)
    return sg_dict


def reverse_complement(sequence):
    trantab = str.maketrans('ACGTacgt', 'TGCAtgca')
    sequence = sequence.translate(trantab)
    return sequence[::-1]


def position(base):
    if base == 'A':
        return 1
    elif base == 'C':
        return 2
    elif base == 'G':
        return 3
    elif base == 'T':
        return 4


def GC(sg):
    seq = sg
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC0_10(sg):
    seq = sg[0:10]
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC10_20(sg):
    seq = sg[10:20]
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC0_5(sg):
    seq = sg[0:5]
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC5_10(sg):
    seq = sg[5:10]
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC10_15(sg):
    seq = sg[10:15]
    return (seq.count('G') + seq.count('C')) / len(seq)


def GC15_20(sg):
    seq = sg[15:20]
    return (seq.count('G') + seq.count('C')) / len(seq)


def gestructure(sg):
    seq = sg.replace('T', 'U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    return mfe_struct


def geenergy(sg):
    seq = sg.replace('T', 'U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    return mfe


def Stem(sg):
    seq = sg.replace('T', 'U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    stem = stru.count('(') * 2
    return stem


def Hairpin(sg):
    seq = sg.replace('T', 'U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    Hairpin = 0
    pattern = re.compile(r'\(\.+\)')
    results = pattern.findall(stru)
    for result in results:
        Hairpin += result.count('.')
    return Hairpin


def Free(sg):
    seq = sg.replace('T', 'U')
    fc = RNA.fold_compound(seq)
    (mfe_struct, mfe) = fc.mfe()
    stru = mfe_struct
    cycle = 0
    pattern = re.compile(r'\(+\.+\)+')
    results = pattern.findall(stru)
    for result in results:
        cycle += len(result)
    free = 20 - cycle
    return free


def CRISPRifeature(sg):
    sg = sg
    for i in range(0, len(sg)):
        value = position(sg[i])
        exec(f'feature_{i + 1} = value')
    gc = GC(sg)
    gc0_10 = GC0_10(sg)
    gc10_20 = GC10_20(sg)
    gc0_5 = GC0_5(sg)
    gc5_10 = GC5_10(sg)
    gc10_15 = GC10_15(sg)
    gc15_20 = GC15_20(sg)
    energy = geenergy(sg)
    stem = Stem(sg)
    hairpin = Hairpin(sg)
    free = Free(sg)
    feature_list = [position(sg[9]), position(sg[19]), position(sg[14]), hairpin, position(sg[5]), gc10_15,
                    position(sg[12]), position(sg[11]), gc10_20, position(sg[3]), stem, position(sg[15]),
                    position(sg[0]), energy, position(sg[18]), gc, position(sg[8]), position(sg[4])]
    return (np.array(feature_list)).reshape(1, 18)


def CRISPRkofeature(sg):
    sg = sg
    gc = GC(sg)
    gc0_10 = GC0_10(sg)
    gc10_20 = GC10_20(sg)
    gc0_5 = GC0_5(sg)
    gc5_10 = GC5_10(sg)
    gc10_15 = GC10_15(sg)
    gc15_20 = GC15_20(sg)
    energy = geenergy(sg)
    stem = Stem(sg)
    hairpin = Hairpin(sg)
    free = Free(sg)
    feature_list = [energy, position(sg[3]), position(sg[19]), position(sg[4]), position(sg[6]), gc10_20, gc0_10,
                    position(sg[13]), position(sg[11]), hairpin, position(sg[1]), gc10_15, gc0_5, position(sg[7]), stem,
                    position(sg[16])]
    return (np.array(feature_list)).reshape(1, 16)


def mismatch(filename):
    mismatch0_dict = {}
    mismatch1_dict = {}
    mismatch2_dict = {}
    mismatch3_dict = {}
    with open(filename, 'r')as offtarget_file:
        content = offtarget_file.read()
    off_list = content.split('\n')
    del off_list[-1]
    for off in off_list:
        sequence = off[0:20]
        if off[-1] == '0':
            mis = off.split('\t')[-3]
            mismatch0_dict.setdefault(sequence, []).append(mis)
            # mismatch0_dict[sequence].append(mis)
        elif off[-1] == '1':
            mis = off.split('\t')[-3]
            mismatch1_dict.setdefault(sequence, []).append(mis)
        elif off[-1] == '2':
            mis = off.split('\t')[-3]
            mismatch2_dict.setdefault(sequence, []).append(mis)
        elif off[-1] == '3':
            mis = off.split('\t')[-3]
            mismatch3_dict.setdefault(sequence, []).append(mis)
    return mismatch0_dict, mismatch1_dict, mismatch2_dict, mismatch3_dict


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


def get_result(res: requests.models.Response,
               csrfmiddlewaretoken: str,
               cookie):
    soup = BeautifulSoup(res.text, 'html.parser')
    try:
        usrid = soup.find('input', {'name': 'usrid'})['value']
    except Exception:
        raise Exception
    params = {
        'csrfmiddlewaretoken': csrfmiddlewaretoken,
        'usrid': usrid
    }
    res = requests.post("http://skl.scau.edu.cn/offtarget/result/",
                        data=params,
                        cookies=cookie)

    soup = BeautifulSoup(res.text, 'html.parser')

    script_tag = soup.find('script')

    js_code = script_tag.text
    find_in_js_code(js_code=js_code, var_declare='var offTargetDict = ')
    find_in_js_code(js_code=js_code, var_declare='var annoteInfoDict = ')
    return find_in_js_code(js_code=js_code, var_declare='var offTargetDict = '), \
           find_in_js_code(js_code=js_code, var_declare='var annoteInfoDict = ')


def find_in_js_code(js_code, var_declare):
    off_target_dict_index = js_code.find(var_declare)
    assert off_target_dict_index != -1
    start_index = off_target_dict_index + len(var_declare)
    end_index = js_code.find(';', start_index)
    off_target_dict_content = js_code[start_index:end_index]
    return json.loads(off_target_dict_content)


class PamMismatchException(Exception):
    pass


# 返回一个字典，键是sgRNA序列（23bp）,值是脱靶点的列表
def process_offtarget_seq(sgrna_seq: str, offtarget_res: str):
    res = sgrna_seq
    cur_idx = 0
    tmp = '0'

    for c in offtarget_res:

        if c.isdigit():
            # print(f"    ~ {tmp} -> {tmp}{c}")
            tmp = f"{tmp}{c}"
        else:
            cur_idx = cur_idx + int(tmp)
            # print(f"    * {cur_idx}{c} {res} => {res[:cur_idx]} {c} {res[cur_idx+1:]}")
            res = res[:cur_idx] + c.lower() + res[cur_idx + 1:]
            tmp = '0'
            cur_idx += 1
    if int(tmp) < 3:
        raise PamMismatchException
    return res


def offtarget_result_process(seqs, result):
    (ot_seq, ot_loc) = result
    # print(f"------- processing ----------")
    # print(f"{ot_seq}\n{ot_loc}\n{len(seqs)}, {len(ot_seq.keys())}, {len(ot_loc.keys())}")
    assert len(seqs) == len(ot_seq.keys()) == len(ot_loc.keys())
    res = dict()
    for i, (origin_seq, origin_pam) in enumerate(seqs):  # 每一个target
        # print('***********')
        index = str(i + 1)
        sgrna_seq = f"{origin_seq}{origin_pam}"
        # print(f"    {index}:{origin_seq} {origin_pam}\n\t{ot_seq[str(index)]}\n\t{ot_loc[str(index)]}\n")
        ot_s = ot_seq[index]
        ot_l = ot_loc[index]
        res[sgrna_seq] = list()
        for k, v in ot_s.items():  # 每一个seq下各个offtarget
            item = v
            try:
                item[3] = process_offtarget_seq(sgrna_seq, item[3])
            except PamMismatchException:
                continue
            item.extend(ot_l[k])
            res[sgrna_seq].append(item)
    return res


def comoff(seqs):
    res = requests.get("http://skl.scau.edu.cn/offtarget")
    if res.status_code == 200:
        soup = BeautifulSoup(res.text, 'html.parser')

        try:
            csrfmiddlewaretoken = soup.find('input', {'name': 'csrfmiddlewaretoken'})['value']
            set_cookie_header = res.headers.get('Set-Cookie', '')
            csrf_token = set_cookie_header.split('csrftoken=')[1].split(';')[0]
        except Exception:
            raise Exception

        payload = {
            'csrfmiddlewaretoken': csrfmiddlewaretoken,
            'pamSelect': 'NGG',
            'otherGuide': 20,
            'genomeSelect': 'GRCh38',
            'seqNum': len(seqs)
        }
        for i, c in enumerate(seqs):
            payload.update({f"seq{i}": c[0], f"pam{i}": c[1]})
        cookie = {'csrftoken': csrf_token}
        jar = RequestsCookieJar()
        for c in set_cookie_header.split(';'):
            key, value = c.split("=", 1)
            jar.set(key, value)

        res = requests.post('http://skl.scau.edu.cn/offtarget/wait/',
                            data=payload,
                            cookies=cookie)

        with open('res.html', 'w', encoding='utf-8') as f:
            f.write(res.text)
        if res.status_code == 200:
            print('首个请求成功！')
            # 可以打印响应内容
            # print(res.text)
            result = get_result(res, csrfmiddlewaretoken, cookie)
            result = offtarget_result_process(seqs, result)
            return result
        else:
            print('首个请求失败:', res)
            return -1
    else:
        print("首个请求失败:", res)
        return -1


# 返回一个字典，键是sgRNA序列（23bp）,值是脱靶点的列表
def getoffsituation(seqs):
    i = 0
    list_length = 20
    off_situation = dict()

    while i < len(seqs):
        result = comoff(seqs[i:i + list_length])
        i = i + list_length
        off_situation.update(result)
    return off_situation  # 返回一个字典，键是sgRNA序列（23bp）,值是脱靶的列表


def binary_search(query, subject):  # query [start, end], subject [[start, end, bio_type, gene], ...]
    query_start, query_end = query
    low = 0
    high = len(subject) - 1
    mid = 0
    while high >= low:
        mid = (low + high) // 2
        # print(mid)
        candidate = subject[mid]
        # print(candidate)
        candidate_start, candidate_end, bio_type = candidate[:3]
        if (
                candidate_start <= query_start and query_end <= candidate_end):  # query total in candidate [candidate_start, query_start, query_end, candidate_end]
            # print(1)
            return candidate  # [start, end, bio_type, gene]
        elif (
                query_start < candidate_start and query_end >= candidate_start and query_end <= candidate_end):  # query_ritht overlap candidate_left
            if bio_type == 'exon':  # if overlap, return exon
                # print(2)
                return candidate
            else:
                # print(3)
                return subject[mid - 1]
        elif (
                candidate_start <= query_start and query_start <= candidate_end and candidate_end < query_end):  # query_left overlap candidate_right
            if bio_type == 'exon':  # if overlap, return exon
                # print(4)
                return candidate
            else:
                # print(5)
                return subject[mid + 1]
        elif (query_start < candidate_start and candidate_end < query_end):  # query_left overlap candidate_right
            if bio_type == 'exon':  # if overlap, return exon
                # print(6)
                return candidate
            else:
                # print(7)
                return subject[mid + 1]
        elif (query_end < candidate_start):  # continue searching on the left
            high = mid - 1
            # print('left')
        elif (query_start > candidate_end):  # continue searching on the right
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
        if (
                candidate_start <= query_start and query_end <= candidate_end):  # query total in candidate [candidate_start, query_start, query_end, candidate_end]
            # print(1)
            return candidate  # [start, end, bio_type, gene]
        elif (
                query_start < candidate_start and query_end >= candidate_start and query_end <= candidate_end):  # query_ritht overlap candidate_left
            if bio_type == 'exon':  # if overlap, return exon
                # print(2)
                return candidate
            else:
                # print(3)
                return subject[mid - 1]
        elif (
                candidate_start <= query_start and query_start <= candidate_end and candidate_end < query_end):  # query_left overlap candidate_right
            if bio_type == 'exon':  # if overlap, return exon
                # print(4)
                return candidate
            else:
                # print(5)
                return subject[mid + 1]
        elif (query_start < candidate_start and candidate_end < query_end):  # query_left overlap candidate_right
            if bio_type == 'exon':  # if overlap, return exon
                # print(6)
                return candidate
            else:
                # print(7)
                return subject[mid + 1]
    return ''  # failure


def getlocation(chrom, strand, startsite, endsite, location_dict):
    location_list = location_dict[chrom]
    location_info = binary_search([startsite, endsite], location_list)
    result_list = []
    if location_info:
        start, end, bio_type, gene = location_info
        # result = [gene, bio_type, '{0}:{1}-{2}'.format(chrom, start, end),strand]
        result = [gene, bio_type]
        result_list.append(result)
    else:
        print('binary search error: {0}:{1}-{2}'.format(chrom, startsite, endsite))
    return result_list[0]


# #返回脱靶得分
# def geoff_score(offsitiation):
#     offscore = {}
#     score_list = []
#     for target in offsitiation:
#         score = 0
#         off_list = offsitiation[target]
#         for off in off_list:
#             mis = count_letters(off[2])
#             if mis == 0:
#                 score = score + 0.6
#             if mis == 1:
#                 score = score + 0.3
#             if mis ==2:
#                 score = score + 0.15
#             if mis >2:
#                 score = score + 0.05
#         offscore[target] = score
#         score_list.append(score)
#         score_nor = normalization(score_list)
#         i = 0
#         for target in offscore:
#             offscore[target] = score_nor[i]
#             i +=1
#     return offscore



#离线版本返回脱靶得分
def geoff_score(offsitiation):
    offscore = {}
    score_list = []
    for target in offsitiation:
        score = 0
        off_list = offsitiation[target]
        for off in off_list:
            mis = count_letters(off)
            if mis == 0:
                score = score + 0.4
            if mis == 1:
                score = score + 0.1
            if mis ==2:
                score = score + 0.07
            if mis >2:
                score = score + 0.05
        offscore[target] = score
        score_list.append(score)
        score_nor = normalization(score_list)
        i = 0
        for target in offscore:
            offscore[target] = score_nor[i]
            i +=1
    return offscore

def normalization(x):
    """"
    归一化到区间{0,1]
    返回副本
    """
    _range = np.max(x) - np.min(x)
    return (x - np.min(x)) / _range


def count_letters(string):
    upper_count = 0
    lower_count = 0
    for i in string:
        if i.islower():
            lower_count += 1
        elif i.isupper():
            upper_count += 1
    return lower_count

def offline_off(filename):
    offtarget_dict = {}
    with open(filename, 'r')as offtarget_file:
        content = offtarget_file.read()
    off_list = content.split('\n')
    del off_list[-1]
    for off in off_list:
        sequence = off[0:20]
        offtarget = off.split('\t')[3]
        if sequence not in offtarget_dict:
            offtarget_dict[sequence] = []
            offtarget_dict[sequence].append(offtarget)
        else:
            offtarget_dict[sequence].append(offtarget)
    return offtarget_dict

if __name__ == '__main__':
    seqs = [
        ('AGTATCTACCAACCTTTGCA', 'AGG'),
        ('AGTATCTACCAACCTATGCA', 'AGG'),
        # ('AGTATCCTACAACCTCCGCA', 'AGG'),
        # ('AGTAGTCTACAACCTTTGCA', 'AGG'),
        # ('AGTATCCTACAACCTTTGCA', 'AGG'),
        # ('AGTATCCTACAACCTTTGCA', 'AGG'),
        # ('AGTATCCTACAACCTTTGCA', 'AGG'),
        # ('AGTATCCTGTAACCTTTGCA', 'AGG'),
        # ('AGTATCCTACAACCTTTGCA', 'AGG'),
        # ('AGTATCCTAAAACCTTTGCA', 'AGG'),
        # ('AGTATCGTTCAACCTTTGCA', 'AGG'),
        # ('AGTATCCTACAAGGTTTGCA', 'AGG'),
        # ('CGTATCTACCGGCCTTTGCA', 'AGG'),
        # ('AGTATCCTACAACCTCCGCA', 'AGG'),
        # ('AGTCCCTTACAACCTTTGCA', 'AGG'),
        # ('GGGCTCCTACAACCTTTGCA', 'AGG'),
        # ('AGGGGGTTACAACCTTTGCA', 'AGG'),
        # ('AGTATCTTTAAACCTTTGCA', 'AGG'),
        # ('ACCGGCCTGTAACCTTTGCA', 'AGG'),
        # ('TTCGATCCTACAACCTTTGC', 'AGG'),
        # ('AGTATCCTACCACTTCCGCA', 'AGG'),
        # ('AGTGTCGTTCAACCTTTGCA', 'AGG'),
    ]
    getoffsituation(seqs)
    # process_offtarget_seq('AGTATCTACCAACCTTTGCAAGG', '4C13T2A0A')
