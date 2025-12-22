import sys, gzip
from itertools import groupby
import regex as re
import time

import argparse
import os


#python3 demutliplex.py fastq.fq.gz .paf index_table output_dir

index_pri_len = 22
index_len = 8

index_find_len = 150

def rev_comp(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

def get_table(index_path, out_path):
    info_dict = {}
    fw_dict = {}
    cnt_dict = {}
    with open(index_path) as f:
        next(f)
        for i in f:
            i = i.strip().split('\t')
            gene = i[0].split('_')[0]
            sid = i[0]
            if gene == 'CD28':
                gene = gene + '_' + i[0][-1]
            if gene not in info_dict:
                info_dict[gene] = [[[i[2][index_pri_len - 10: index_pri_len], i[2][index_pri_len + index_len: index_pri_len + index_len + 10]],  [rev_comp(i[4][index_pri_len + index_len: index_pri_len + index_len+10]), rev_comp(i[4][index_pri_len - 10: index_pri_len])]],{}]
            #info_dict[gene][sid] = [i[2][index_pri_len+index_len:], i[3], i[4][index_pri_len+index_len:], i[5]]
            info_dict[gene][1][sid] = [i[3],  rev_comp(i[5])]
            fw_dict[sid] = open(f'{out_path}/{sid}.fastq', 'w')
            cnt_dict[sid] = 0
    fw_dict['others']  = open(f'{out_path}/others.fastq', 'w')
    cnt_dict['others'] = 0

    return info_dict, fw_dict, cnt_dict

def fastq_iter(path):
    with gzip.open(path) as f:
        while True:
            h = f.readline().decode('utf-8')
            if not h: break
            if not h.startswith('@'):
                raise RuntimeError("FASTQ 형식이 아님")
            rid = h[1:].strip().split()[0]
            seq = f.readline().strip().decode('utf-8')
            f.readline()  
            qual = f.readline().strip().decode('utf-8')
            yield rid, seq, qual

def paf_group_iter(path):
    def key(line): return line.split('\t',1)[0]
    with open(path) as f:
        for qname, grp in groupby(f, key=key):
            yield qname, list(grp)



def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 9999  
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def match_sample(db: dict, query_f: str, query_b: str):
    for name, (f, b) in db.items():
        if f == query_f and b == query_b:
            return name
    candidates = []
    for name, (f, b) in db.items():
        if hamming(f, query_f) <= 1 and hamming(b, query_b) <= 1:
            candidates.append(name)
    return candidates[0] if len(candidates) == 1 else None

def compare_pri(seq, pri_seq, front=True):
    for err_n in range(1):

        pat = re.compile(rf"(?b)({pri_seq[0]}){{e<={err_n}}}([ACGTN]{{8}})({pri_seq[1]}){{e<={err_n}}}", flags=re.IGNORECASE)
        res = pat.finditer(seq, overlapped=True)
        res= list(res)

        if res != []:
            if front:
                r = res[0]
                idx = r.group(2)
                return idx,  r.end()
            else:
                r = res[-1]
                idx = r.group(2)
                return idx,  r.start()
    return None, None


def find_index(reads, info, info_dict):
    info = info.strip().split('\t')
    ref_name = info[5]
    strand = info[4]
    st = int(info[2]) - index_find_len + 20 
    ed = int(info[3]) + index_find_len - 20
    if st < 0 : st = 0
    if ed > int(info[1]): ed = int(info[1])
    s = reads[1][st: ed]
    q = reads[2][st: ed]
    if len(s) < 300:
        return False
    if strand == '-':
        s = rev_comp(s)
        q = q[::-1]

    front_idx, front_end = compare_pri(s[:index_find_len], info_dict[ref_name][0][0])
    back_idx, back_start = compare_pri(s[-index_find_len:], info_dict[ref_name][0][1], front=False)

    #print(reads[0])
    #print(s[:index_find_len])
    #print(s[-index_find_len:])
    #print(front_idx, back_idx)
    
    if front_idx == None or back_idx == None:
        return False
    back_start = -index_find_len + back_start - 10
    front_end -= 10
    sample_name = match_sample(info_dict[ref_name][1], front_idx, back_idx)

    #print(sample_name)
    #input()

    if sample_name == None:
        return 'others', reads[0], s[front_end: back_start], q[front_end: back_start]

    return sample_name, reads[0], s[front_end: back_start], q[front_end: back_start]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Demultiplex FASTQ reads using PAF alignment and index table"
    )

    # positional arguments
    parser.add_argument(
        "fastq",
        help="Input gzipped FASTQ file"
    )
    parser.add_argument(
        "paf",
        help="PAF alignment file (from minimap2)"
    )
    parser.add_argument(
        "index_table",
        help="Index table mapping barcodes to sample IDs"
    )
    parser.add_argument(
        "output_dir",
        help="Output directory"
    )

    # optional arguments
    parser.add_argument(
        "--index_pri_len",
        type=int,
        default=22,
        help="Length of primer sequence flanking the index (default: 22)"
    )
    parser.add_argument(
        "--index_len",
        type=int,
        default=8,
        help="Length of barcode/index sequence (default: 8)"
    )
    parser.add_argument(
        "--index_find_len",
        type=int,
        default=150,
        help="Search window size from alignment boundary (default: 150 bp)"
    )

    return parser.parse_args()


def main():

    global index_pri_len, index_len, index_find_len

    args = parse_args()

    # sanity checks
    if not os.path.exists(args.fastq):
        raise FileNotFoundError(f"FASTQ not found: {args.fastq}")
    if not os.path.exists(args.paf):
        raise FileNotFoundError(f"PAF not found: {args.paf}")
    if not os.path.exists(args.index_table):
        raise FileNotFoundError(f"Index table not found: {args.index_table}")

    os.makedirs(args.output_dir, exist_ok=True)

    print("=== Demultiplexing parameters ===")
    print(f"FASTQ           : {args.fastq}")
    print(f"PAF             : {args.paf}")
    print(f"Index table     : {args.index_table}")
    print(f"Output dir      : {args.output_dir}")
    print(f"Primer length   : {args.index_pri_len}")
    print(f"Index length    : {args.index_len}")
    print(f"Search window   : ±{args.index_find_len} bp")
    print("================================")

    info_dict, fw_dict, cnt_dict = get_table(args.index_table, args.output_dir)

    fq = fastq_iter(args.fastq)
    pafg = paf_group_iter(args.paf)

    index_pri_len = args.index_pri_len
    index_len = args.index_len
    index_find_len = args.index_find_len

    fq_reads = next(fq)
    paf_key, paf_lines = next(pafg)

    line_n = 0
    t = time.time()

    while True:

        if paf_key == None or fq_reads[0] == None:
            break
        if paf_key != fq_reads[0]:
            print('unmatched! paf_key, fq_reads[0]')
            while paf_key != fq_reads[0]:
                fq_reads = next(fq, None)
                if fq_reads == None:
                    break
            print('Found')
            if fq_reads == None:
                break
        line_n += 1

        paf_n = 0
        for paf_res in paf_lines:
            paf_n += 1
            sample_info = find_index(fq_reads, paf_res, info_dict)
            if sample_info == False:
                continue
            fw_dict[sample_info[0]].write(f'@{sample_info[1]}_{paf_n}\n{sample_info[2]}\n+\n{sample_info[3]}\n')
            cnt_dict[sample_info[0]] += 1
        if line_n % 5000 == 0:
            print(f'Demultiplexing {line_n} reads... {round(time.time() - t,2)}         \r', end ='')
            t = time.time()
        

        fq_reads = next(fq, None)
        paf_key, paf_lines = next(pafg, (None, None))
        if paf_key == None:
            break

    for x,y in fw_dict.items():
        y.close()

    fw = open('summary.txt', 'w')
    for x, y in cnt_dict.items():
        fw.write(f'{x}\t{y}\n')
        print(f'{x}\t{y}')
    fw.close()

if __name__ == "__main__":
    main()


