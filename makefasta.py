import re
import argparse

def get_args():
  parser = argparse.ArgumentParser(description="make fasta from list")
  parser.add_argument("-list")
  return parser.parse_args()

args = get_args()
list = args.list

out = open("output_fasta.fa","a+")
with open(list, "r+") as fh:
    for line in fh:
        out.write(">")
        out.write("\n")
        out.write(line)
