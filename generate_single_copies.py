import re
import argparse
import numpy
import random
import math

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-array")
  parser.add_argument("-assembly")
  parser.add_argument("-sc")
  parser.add_argument("-k")
  parser.add_argument("-asp")
  return parser.parse_args()

args = get_args()
array = args.array
assembly = args.assembly
single_copy = args.sc
array_specific_set = args.asp
k = int(args.k)

def gc_content(sequence):
    AT = 0
    GC = 0
    No_Ns = True
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
            No_Ns = False
            break
    if No_Ns == True:
        total = GC+AT
        ratio = GC/total
        return(round(ratio,2))
    else:
        return("Null")


with open(array_specific_set, "r+") as set:
    kmer_dic = {}
    for line in set:
        if line[0] == "@":
            regex = re.search(r'@([\S]+)',line)
            Array_GC = float(regex.group(1))
        elif line[0] == ">":
            regex = re.search(r'\>([\S]+)',line)
            chrID = regex.group(1)
        else:
            regex = re.search(r'([\S]+)\t([\S]+)',line)
            kmer = regex.group(1)
            count = int(regex.group(2))
            kmer_dic[kmer] = count


#sum array copies to use in the single copy matching
copy_number = 0
for key in kmer_dic:
    copy_number += kmer_dic[key]

print("Scanning rest of chromosome for 10kb matched G/C regions")
valid_window_list = []
with open(assembly, "r+") as ref:
    while True:
        seqname = ref.readline().strip()
        if seqname == "":
            break
        seq = ref.readline().strip()
        seqID = re.search(r'\>([\S]+)\:[\S]+',seqname)
        chrmID = seqID.group(1)
        if chrmID == chrID:
            for num in range(int(math.floor(len(seq)/10000))):
                window = seq[(0+10000*num):(10000-1+10000*num)]
                GC_ratio = gc_content(window)
                if GC_ratio != "Null":
                    if abs(Array_GC - GC_ratio) < 0.03:
                        valid_window_list.append([0+10000*num,10000-1+10000*num])


single_copy_dic = {}
for n in range(101):
    ratio = n/100
    single_copy_dic[ratio] = []


print("Selecting single copy kmers from G/C compatible genomic regions")
print(len(kmer_dic))
print(copy_number)

incorrect_single_copy_kmers = []
with open(single_copy) as scfa:
    i = 0
    while True:
        chrmname = scfa.readline().strip()
        if chrmname == "":
            break
        seq = scfa.readline().strip()
        seqID = re.search(r'\>([\S]+)\:([\S]+)\-([\S]+)',chrmname)
        chrmID = seqID.group(1)
        start_pos = int(seqID.group(2))
        end_pos = int(seqID.group(3))
        if chrmID == chrID:
            for index in range(len(valid_window_list)):
                start = valid_window_list[index][0]
                end = valid_window_list[index][1]
                if abs(start_pos-start) < abs(start-end) & abs(start_pos - end) < abs(start-end):
                    if abs(end_pos-start) < abs(start-end) & abs(end_pos - end) < abs(start-end):
                        gc_ratio = gc_content(seq)
                        if seq not in single_copy_dic[gc_ratio]:
                            single_copy_dic[gc_ratio].append(seq)
                            i += 1
                            if i % 50000 == 0:
                                print("working")
                            break
        #if i > (copy_number*100):
        if i > (copy_number*5):
            print("should be stopping soon")
            break


# total_len = 0
# for index in single_copy_dic:
#     print(single_copy_dic[index])
#     a <- list(single_copy_dic[index])
#     b <- len(list(single_copy_dic[index]))
#     total_len += b
#     if len(set(list(single_copy_dic[index]))) != len(list(single_copy_dic[index])):
#         print("Oh No! The single copy kmers are duplicated")


print("Matching GC ratio between array-specific and single-copy kmers")
single_copy_set = {}
key_list = list(single_copy_dic.keys())
val_list = list(single_copy_dic.values())

#print(len(numpy.unique(key_list)))
#print(len(kmer_dic.keys()))

# print(total_len)
print(copy_number)
total_copies = 0
for item in kmer_dic:
    gc_ratio = gc_content(item)
    copy_number = kmer_dic[item]
    total_copies += copy_number
    #print(copy_number)
    for num in range(copy_number):
        if single_copy_dic[gc_ratio] != []:
            kmer = random.choice(single_copy_dic[gc_ratio])
            #print(kmer)
            if kmer in single_copy_set:
                print("What!?")
            single_copy_set[kmer] = 0
            #del single_copy_dic[gc_ratio][kmer]
            #print(single_copy_dic[gc_ratio].index(kmer))
            #print(single_copy_dic[gc_ratio][single_copy_dic[gc_ratio].index(kmer)])
            #print(single_copy_dic[gc_ratio][single_copy_dic[gc_ratio].index(kmer)])
            del single_copy_dic[gc_ratio][single_copy_dic[gc_ratio].index(kmer)]
        else:
            min_dist = 1.1
            for value in single_copy_dic.keys():
                dist = abs(gc_ratio - value)
                if dist < min_dist:
                    if single_copy_dic[value] != []:
                        min_dist = dist
                        min_dist_value = value
            kmer = random.choice(single_copy_dic[min_dist_value])
            #print(kmer)
            if kmer in single_copy_set:
                print("What!?")
            single_copy_set[kmer] = 0
            #print(single_copy_dic[gc_ratio][single_copy_dic[gc_ratio].index(kmer)])
            del single_copy_dic[min_dist_value][single_copy_dic[min_dist_value].index(kmer)]
#print(single_copy_set)


print(len(single_copy_set))
print(total_copies)

set = open("single_copy_set.tsv","a+")
for entry in single_copy_set:
    set.write(entry)
    set.write("\t")
    set.write(str(single_copy_set[entry]))
    set.write("\n")
