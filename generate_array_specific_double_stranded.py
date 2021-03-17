import re
import argparse
import numpy
import random
import math

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-array")
  parser.add_argument("-assembly")
  parser.add_argument("-k")
  # parser.add_argument("-copycutoff")
  # parser.add_argument("-less")
  return parser.parse_args()

args = get_args()
array = args.array
assembly = args.assembly
k = int(args.k)
# moreless = args.less
# cutoff = args.copycutoff

def gc_content(sequence):
    AT = 0
    GC = 0
    for character in sequence:
        if character == "A":
            AT += 1
        elif character == "T":
            AT += 1
        elif character == "G":
            GC += 1
        elif character == "C":
            GC += 1
        else:
            print("Non-uppercase ATGC character detected")
            break
    total = GC+AT
    ratio = GC/total
    return(round(ratio,2))

def reverse_complement(sequence):
    new_seq = ""
    for character in sequence:
        if character == "A":
            new_seq += "T"
        elif character == "T":
            new_seq += "A"
        elif character == "G":
            new_seq += "C"
        elif character == "C":
            new_seq += "G"
        else:
            print("Non-uppercase ATGC character detected")
            break
    reversed_seq = new_seq[::-1]
    return(new_seq)


print("Generating array kmers")
with open(array, "r+") as fh:
    kmer_dic = {}
    seqname = fh.readline().strip()
    seq = fh.readline().strip()
    seqID = re.search(r'\>([\S]+)\:[\S]+',seqname)
    chrID = seqID.group(1)
    Array_GC = gc_content(seq)
    for n in range(len(seq)-k+1):
        kmer = seq[n:n+k]
        if kmer in kmer_dic:
            kmer_dic[kmer] += 1
        else:
            kmer_dic[kmer] = 1

# below_cutoff = []
# if moreless == "more":
#     for item in kmer_dic:
#         if int(kmer_dic[item]) < int(cutoff):
#             below_cutoff.append(item)
# elif moreless == "less":
#     for item in kmer_dic:
#         if int(kmer_dic[item]) > int(cutoff):
#             below_cutoff.append(item)
# for kmer in below_cutoff:
#     del kmer_dic[kmer]

print("Identifying array-specific kmers")
with open(assembly, "r+") as ref:
    for line in ref:
        seq = line.strip()
        if seq[0] == ">":
            continue
        kmer_list = []
        for n in range(len(seq)-k+1):
            kmer = seq[n:n+k]
            kmer_revcomp = reverse_complement(kmer)
            if kmer in kmer_dic:
                #print(kmer)
                del kmer_dic[kmer]
            elif kmer_revcomp in kmer_dic:
                del kmer_dic[kmer_revcomp]



# graph = open("kmer_frequency_x_array.tsv","a+")
# for entry in kmer_dic:
#     graph.write(entry+"\t"+str(kmer_dic[entry]))
#     graph.write("\n")

set = open("array_specific_set.tsv","a+")
set.write("@"+str(Array_GC))
set.write("\n")
set.write(">"+str(chrID))
set.write("\n")
for entry in kmer_dic:
    set.write(entry)
    set.write("\t")
    set.write(str(kmer_dic[entry]))
    set.write("\n")
