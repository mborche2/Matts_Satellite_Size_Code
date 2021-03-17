import re
import argparse
import numpy
import random
import math

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-wgs1")
  parser.add_argument("-wgs2")
  parser.add_argument("-k")
  parser.add_argument("-asp")
  parser.add_argument("-scs")
  parser.add_argument("-ras")
  return parser.parse_args()

args = get_args()
reads1 = args.wgs1
reads2 = args.wgs2
array_specific = args.asp
single_copy = args.scs
k = int(args.k)
reference_array_size = int(args.ras)

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

with open(array_specific, "r+") as set:
    kmer_dic = {}
    for line in set:
        if line[0] == "@":
            regex = re.search(r'@([\S]+)',line)
            Array_GC = regex.group(1)
        elif line[0] == ">":
            continue
        else:
            regex = re.search(r'([\S]+)\t([\S]+)',line)
            kmer = regex.group(1)
            count = int(regex.group(2))
            kmer_dic[kmer] = count

with open(single_copy, "r+") as set:
    single_copy_set = {}
    for line in set:
        regex = re.search(r'([\S]+)\t([\S]+)',line)
        kmer = regex.group(1)
        count = int(regex.group(2))
        single_copy_set[kmer] = count


kmer_cov_dic = {}
for key in kmer_dic:
    kmer_cov_dic[key] = 0

with open(reads1, "r+") as fil:
    i = 0
    for line in fil:
        i += 1
        if i%4 == 2:
            read = line.strip()
        if i%4 ==0:
            use = True
            line = line.strip()
            for char in line:
                score = ord(char) - 33
                if score < 21:
                    use = False
            if use == True:
                kmer_list = []
                for n in range(len(read)-k+1):
                    kmer = read[n:n+k]
                    kmer_comp = reverse_complement(kmer)
                    if kmer in kmer_cov_dic:
                        kmer_cov_dic[kmer] += 1
                    elif kmer_comp in kmer_cov_dic:
                        kmer_cov_dic[kmer_comp] += 1
                    elif kmer in single_copy_set:
                        single_copy_set[kmer] += 1
                    elif kmer_comp in single_copy_set:
                        single_copy_set[kmer_comp] += 1

with open(reads2, "r+") as fil:
    i = 0
    for line in fil:
        i += 1
        if i%4 == 2:
            read = line.strip()
        if i%4 ==0:
            use = True
            line = line.strip()
            for char in line:
                score = ord(char) - 33
                if score < 21:
                    use = False
            if use == True:
                kmer_list = []
                for n in range(len(read)-k+1):
                    kmer = read[n:n+k]
                    kmer_comp = reverse_complement(kmer)
                    if kmer in kmer_cov_dic:
                        kmer_cov_dic[kmer] += 1
                    elif kmer_comp in kmer_cov_dic:
                        kmer_cov_dic[kmer_comp] += 1
                    elif kmer in single_copy_set:
                        single_copy_set[kmer] += 1
                    elif kmer_comp in single_copy_set:
                        single_copy_set[kmer_comp] += 1


#calculation equation
#(chm13 array size) * (summed coverage of all array-specific k-mers)/(total copy number of array-specific k-mers in chm13))/(summed coverage over some reference set of single-copy k-mers that are AT-matched and on the same chr))

#summed coveraged of array specific kmers
coverage_sum = 0
for key in kmer_cov_dic:
    coverage_sum += kmer_cov_dic[key]
print("The array specific coverage sum:")
print(coverage_sum)

#summed coverage of AT matched single copy kmers
sc_coverage_sum = 0
for key in single_copy_set:
    sc_coverage_sum += single_copy_set[key]
print("The single copy coverage sum:")
print(sc_coverage_sum)

#total copy number array specific kmers
copy_number = 0
for key in kmer_dic:
    copy_number += kmer_dic[key]
print("The total number of kmers in the array-specific set:")
print(copy_number)

kmer_calculated_array_length = reference_array_size*(coverage_sum)/sc_coverage_sum
#kmer_calculated_array_length = 3106918*(coverage_sum/copy_number)/sc_coverage_sum
print("The calculated array length!:")
print(kmer_calculated_array_length)

output_file = open("array_size.txt","a+")
output_file.write("The array specific coverage sum:")
output_file.write(str(coverage_sum))
output_file.write("\n")
output_file.write("The single copy coverage sum:")
output_file.write(str(sc_coverage_sum))
output_file.write("\n")
output_file.write("The total number of kmers in the array-specific set:")
output_file.write(str(copy_number))
output_file.write("\n")
output_file.write("The calculated array length!:")
output_file.write(str(kmer_calculated_array_length))
output_file.write("\n")
