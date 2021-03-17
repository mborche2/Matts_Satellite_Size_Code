import re
import argparse
import numpy
import random
import math
import multiprocessing as mp

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

def loop_kmers(fil):
    with open(fil) as scfa:
        i = 0
        n=0
        first_line = True
        while True:
            chrmname = scfa.readline().strip()
            # n += 1
            # if n > 20000:
            #     print("stopping no sc detected for test")
            #     break
            if first_line == True:
                if chrmname == "":
                    chrmname = scfa.readline().strip()
                    if chrmname[0] != ">":
                        chrmname = scfa.readline().strip()
                        if chrmname[0] != ">":
                            chrmname = scfa.readline().strip()
                            print(chrmname)
                            first_line = False
                        else:
                            first_line = False
                    else:
                        first_line = False
                elif chrmname[0] != ">":
                    chrmname = scfa.readline().strip()
                    if chrmname[0] != ">":
                        chrmname = scfa.readline().strip()
                        print(chrmname)
                        first_line = False
                    else:
                        first_line = False
                else:
                    first_line = False
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
            if i > (copy_number*1):
                print("should be stopping soon")
                break
    return(single_copy_dic)

def start_process():
    print('Starting', mp.current_process().name)

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
        if seqID.group(1) == False:
            print(seqname)
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

# with open(single_copy) as scfa:
#     length = 0
#     for line in scfa:
#         length += 1
#     scfa.seek(0)
#     scfa.close()


incorrect_single_copy_kmers = []
chromosome = chrID.strip("chr")
file1 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/single_copy_loci/single_copy_by_chrm_shuffled/" + chromosome + "_shufaa"
file2 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/single_copy_loci/single_copy_by_chrm_shuffled/" + chromosome + "_shufab"
file3 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/single_copy_loci/single_copy_by_chrm_shuffled/" + chromosome + "_shufac"
file4 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/single_copy_loci/single_copy_by_chrm_shuffled/" + chromosome + "_shufad"
file5 = "/n/core/bigDataAI/Genomics/Gerton/jennifer_gerton/jeg10/single_copy_loci/single_copy_by_chrm_shuffled/" + chromosome + "_shufae"
file_list = [file1,file2,file3,file4,file5]

# pool = mp.Pool(5)
# pool.map(loop_kmers, [fil for fil in file_list])
# pool.close()

if __name__ == '__main__':
    pool = mp.Pool(processes= 5,initializer=start_process)
    pool_outputs1 = pool.map(loop_kmers, [fh for fh in file_list])
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

print("The number of pool outputs is:")
print(len(pool_outputs1))

single_copy_set_dic_list=numpy.empty(5,dtype='object')
for i in range(5):
    single_copy_set_dic_list[i] = pool_outputs1[i]


single_copy_dic={}

for n in range(101):
    ratio = n/100
    temp_list=[]
    for i in range(5):
        temp_list = temp_list + single_copy_set_dic_list[i][ratio]
    single_copy_dic[ratio] = temp_list


#check if any kmers occur multiple times
# for index in single_copy_dic:
#     if len(set(single_copy_dic)) != len(single_copy_dic):
#         print("Oh No! The single copy kmers are duplicated")
#         print(exit)

print("Matching GC ratio between array-specific and single-copy kmers")
single_copy_set = {}
key_list = list(single_copy_dic.keys())
val_list = list(single_copy_dic.values())

#print(len(numpy.unique(key_list)))
#print(len(kmer_dic.keys()))


for item in kmer_dic:
    gc_ratio = gc_content(item)
    copy_number = kmer_dic[item]
    #print(copy_number)
    for num in range(copy_number):
        if single_copy_dic[gc_ratio] != []:
            kmer = random.choice(single_copy_dic[gc_ratio])
            if kmer in single_copy_set:
                print("WTF")
            single_copy_set[kmer] = 0
            #del single_copy_dic[gc_ratio][kmer]
            #print(single_copy_dic[gc_ratio].index(kmer))
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
            if kmer in single_copy_set:
                print("WTF")
            single_copy_set[kmer] = 0
            del single_copy_dic[min_dist_value][single_copy_dic[min_dist_value].index(kmer)]
#print(single_copy_set)


print(len(single_copy_set))

set = open("single_copy_set.tsv","a+")
for entry in single_copy_set:
    set.write(entry)
    set.write("\t")
    set.write(str(single_copy_set[entry]))
    set.write("\n")
