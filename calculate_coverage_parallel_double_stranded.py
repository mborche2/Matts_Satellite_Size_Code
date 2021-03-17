import re
import argparse
import numpy
import random
import math
import multiprocessing as mp
from multiprocessing import Process, Manager
import operator
from collections import Counter

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

def loop_reads(fh):
    with open(fh, "r+") as fil:
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
    return([kmer_cov_dic,single_copy_set])

def start_process():
    print('Starting', mp.current_process().name)

def combine_dicts(a, b, c, d, e, f, g, h, i, j, op=operator.add):
    return dict(list(a.items()) + list(b.items()) + list(c.items()) + list(d.items()) + list(e.items()) + list(f.items()) + list(g.items()) + list(h.items()) + list(i.items()) + list(j.items()) +
        [(k, op(a[k], b[k], c[k], d[k], e[k], f[k], g[k], h[k], i[k], j[k])) for k in set(j) & set(i) & set(h) & set(g) & set(f) & set(e) & set(d) & set(c) & set(b) & set(a)])


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
    return(reversed_seq)


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

kmer_cov_dic = Counter({})
single_copy_set = Counter({})
for key in kmer_dic:
    kmer_cov_dic[key] = 0
with open(single_copy, "r+") as set:
    for line in set:
        regex = re.search(r'([\S]+)\t([\S]+)',line)
        kmer = regex.group(1)
        count = int(regex.group(2))
        single_copy_set[kmer] = count


file1 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "aa"
file2 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" +  "ab"
file3 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ac"
file4 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ad"
file5 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ae"
file_list = [file1,file2,file3,file4,file5]

file1 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "aa"
file2 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ab"
file3 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ac"
file4 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ad"
file5 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ae"
file_list2 = [file1,file2,file3,file4,file5]


if __name__ == '__main__':
    pool = mp.Pool(processes= 5,initializer=start_process)
    pool_outputs1 = pool.map(loop_reads, [fh for fh in file_list])
    pool.close()
    pool.join()
    # with Pool(5) as p:
    #     for i in range(5):
    #         print(p.map(loop_reads, file_list[i]))


# if __name__ == '__main__':
#     q2 = mp.Queue()
#     for i in range(5):
#         p2 = mp.Process(target=loop_reads, args=(file_list2[i], q2,))
#         p2.start()

kmer_cov_dic_list=numpy.empty(5,dtype='object')
single_copy_set_dic_list=numpy.empty(5,dtype='object')
for i in range(5):
    kmer_cov_dic_list[i], single_copy_set_dic_list[i] = pool_outputs1[i]

if __name__ == '__main__':
    pool = mp.Pool(processes= 5,initializer=start_process)
    pool_outputs2 = pool.map(loop_reads, [fh for fh in file_list2])
    pool.close()
    pool.join()

kmer_cov_dic_list2=numpy.empty(5,dtype='object')
single_copy_set_dic_list2=numpy.empty(5,dtype='object')
for i in range(5):
    kmer_cov_dic_list2[i], single_copy_set_dic_list2[i] = pool_outputs2[i]

sum1= 0
sum2= 0
sum3=0
sum4=0
sum5=0

sum6= 0
sum7= 0
sum8=0
sum9=0
sum10=0

#test if values from
for entry in kmer_cov_dic_list[0]:
    sum1 += kmer_cov_dic_list[0][entry]

for entry in kmer_cov_dic_list[1]:
    sum2 += kmer_cov_dic_list[1][entry]

for entry in kmer_cov_dic_list[2]:
    sum3 += kmer_cov_dic_list[2][entry]

for entry in kmer_cov_dic_list[3]:
    sum4 += kmer_cov_dic_list[3][entry]

for entry in kmer_cov_dic_list[4]:
    sum5 += kmer_cov_dic_list[4][entry]

for entry in single_copy_set_dic_list[0]:
    sum6 += single_copy_set_dic_list[0][entry]

for entry in single_copy_set_dic_list[1]:
    sum7 += single_copy_set_dic_list[1][entry]

for entry in single_copy_set_dic_list[2]:
    sum8 += single_copy_set_dic_list[2][entry]

for entry in single_copy_set_dic_list[3]:
    sum9 += single_copy_set_dic_list[3][entry]

for entry in single_copy_set_dic_list[4]:
    sum10 += single_copy_set_dic_list[4][entry]

# for entry in kmer_cov_dic_list2[0]:
#     sum3 += kmer_cov_dic_list2[0][entry]
#
# for entry in kmer_cov_dic_list2[1]:
#     sum4 += kmer_cov_dic_list2[1][entry]

print("sum of process 1, read 1")
print(sum1)
print("sum of process 2, read 1")
print(sum2)
print("sum of process 3, read 1")
print(sum3)
print("sum of process 4, read 1")
print(sum4)
print("sum of process 5, read 1")
print(sum5)

print("sc sum of process 1, read 1")
print(sum6)
print("sc sum of process 2, read 1")
print(sum7)
print("sc sum of process 3, read 1")
print(sum8)
print("sc sum of process 4, read 1")
print(sum9)
print("sc sum of process 5, read 1")
print(sum10)
# print("sum of process 1, read 2")
# print(sum3)
# print("sum of process 2, read 2")
# print(sum4)

merged_kmer_cov_dic=kmer_cov_dic_list[0]+kmer_cov_dic_list[1]+kmer_cov_dic_list[2]+kmer_cov_dic_list[3]+kmer_cov_dic_list[4]+kmer_cov_dic_list2[0]+kmer_cov_dic_list2[1]+kmer_cov_dic_list2[2]+kmer_cov_dic_list2[3]+kmer_cov_dic_list2[4]
merged_single_copy_set=single_copy_set_dic_list[0]+single_copy_set_dic_list[1]+single_copy_set_dic_list[2]+single_copy_set_dic_list[3]+single_copy_set_dic_list[4]+single_copy_set_dic_list2[0]+single_copy_set_dic_list2[1]+single_copy_set_dic_list2[2]+single_copy_set_dic_list2[3]+single_copy_set_dic_list2[4]
#merged_kmer_cov_dic = pool_outputs1[0][0] + pool_outputs1[1][0] + pool_outputs1[2][0] + pool_outputs1[3][0] + pool_outputs1[4][0] + pool_outputs2[0][0] + pool_outputs2[1][0] + pool_outputs2[2][0] + pool_outputs2[3][0] + pool_outputs2[4][0]
#merged_single_copy_set = pool_outputs1[0][1] + pool_outputs1[1][1] + pool_outputs1[2][1] + pool_outputs1[3][1] + pool_outputs1[4][1] + pool_outputs2[0][1] + pool_outputs2[1][1] + pool_outputs2[2][1] + pool_outputs2[3][1] + pool_outputs2[4][1]


# merged_kmer_cov_dic = combine_dicts(pool_outputs1[0][0],pool_outputs1[1][0],pool_outputs1[2][0],pool_outputs1[3][0],pool_outputs1[4][0],pool_outputs2[0][0],pool_outputs2[1][0],pool_outputs2[2][0],pool_outputs2[3][0],pool_outputs2[4][0],operator.add)
# merged_single_copy_set = combine_dicts(pool_outputs1[0][1],pool_outputs1[1][1],pool_outputs1[2][1],pool_outputs1[3][1],pool_outputs1[4][1],pool_outputs2[0][1],pool_outputs2[1][1],pool_outputs2[2][1],pool_outputs2[3][1],pool_outputs2[4][1],operator.add)

# file1 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "aa"
# file2 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ab"
# file3 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ac"
# file4 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ad"
# file5 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split1" + "ae"
# file_list = [file1,file2,file3,file4,file5]
#
# pool = mp.Pool(5)
# pool.map(loop_reads, [fh for fh in file_list])
# pool.close()
#
# file1 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "aa"
# file2 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ab"
# file3 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ac"
# file4 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ad"
# file5 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/WGS/Chm13/" + "temp_split2" + "ae"
# file_list = [file1,file2,file3,file4,file5]
#
# pool = mp.Pool(5)
# pool.map(loop_reads, [fh for fh in file_list])
# pool.close()


#calculation equation
#(chm13 array size) * (summed coverage of all array-specific k-mers)/(total copy number of array-specific k-mers in chm13))/(summed coverage over some reference set of single-copy k-mers that are AT-matched and on the same chr))

#summed coveraged of array specific kmers

#coverage_sum1 = 0
#fake_kmer_cov_dic = pool_outputs1[1][0]
#for key in fake_kmer_cov_dic:
#    coverage_sum1 += fake_kmer_cov_dic[key]
#print(coverage_sum1)
#
#coverage_sum2 = 0
#fake_scs_dic = pool_outputs1[1][0]
#for key in fake_kmer_cov_dic:
#    coverage_sum2 += fake_scs_dic[key]
#print(coverage_sum2)

coverage_sum = 0
for key in merged_kmer_cov_dic:
    coverage_sum += merged_kmer_cov_dic[key]
print("The array specific coverage sum:")
print(coverage_sum)

#summed coverage of AT matched single copy kmers
sc_coverage_sum = 0
for key in merged_single_copy_set:
    sc_coverage_sum += merged_single_copy_set[key]
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
