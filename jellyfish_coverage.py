
import re
import argparse
import numpy
import random
import math
from collections import Counter

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-asp1")
  parser.add_argument("-asp2")
  parser.add_argument("-k")
  parser.add_argument("-scs1")
  parser.add_argument("-scs2")
  parser.add_argument("-ras")
  return parser.parse_args()

args = get_args()
asp_counts1 = args.asp1
asp_counts2 = args.asp2
scs_counts1 = args.scs1
scs_counts2 = args.scs2
k = int(args.k)
reference_array_size = int(args.ras)

with open(asp_counts1, "r+") as set:
    asp1_counter = Counter({})
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmercounts = countline.strip(">")
        kmer = set.readline().strip()
        asp1_counter[kmer] = int(kmercounts)

with open(asp_counts2, "r+") as set:
    asp2_counter = Counter({})
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmercounts = countline.strip(">")
        kmer = set.readline().strip()
        asp2_counter[kmer] = int(kmercounts)

with open(scs_counts1, "r+") as set:
    scs1_counter = Counter({})
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmercounts = countline.strip(">")
        kmer = set.readline().strip()
        scs1_counter[kmer] = int(kmercounts)

with open(scs_counts2, "r+") as set:
    scs2_counter = Counter({})
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmercounts = countline.strip(">")
        kmer = set.readline().strip()
        scs2_counter[kmer] = int(kmercounts)

kmer_cov_dic = asp1_counter + asp2_counter
single_copy_set = scs1_counter + scs2_counter

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
# copy_number = 0
# for key in kmer_dic:
#     copy_number += kmer_dic[key]
# print("The total number of kmers in the array-specific set:")
# print(copy_number)

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
# output_file.write("The total number of kmers in the array-specific set:")
# output_file.write(str(copy_number))
# output_file.write("\n")
output_file.write("The calculated array length!:")
output_file.write(str(kmer_calculated_array_length))
output_file.write("\n")
