import re
import argparse
import numpy
import random
import math
from collections import Counter

def get_args():
  parser = argparse.ArgumentParser(description="Find kmers that uniquely define a satellite array in a genome")
  parser.add_argument("-array")
  parser.add_argument("-assembly")
  return parser.parse_args()

args = get_args()
arr = args.array
assm = args.assembly


with open(arr, "r+") as set:
    arr_counter = {}
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmercounts = countline.strip(">")
        kmer = set.readline().strip()
        arr_counter[kmer] = int(kmercounts)

with open(assm, "r+") as set:
    while True:
        countline = set.readline().strip()
        if countline == "":
            break
        kmer = set.readline().strip()
        if kmer in arr_counter.keys():
            del arr_counter[kmer]

output_file = open("jellyfish_array_specific_set.txt","a+")
for item in arr_counter:
    output_file.write(item)
    output_file.write("\t")
    output_file.write(str(arr_counter[item]))
    output_file.write("\n")
