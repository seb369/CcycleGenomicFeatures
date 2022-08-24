#!/home/sam/anaconda3/bin/python

# Writen by Roland Wilhelm and modified by Sam Barnett
import sys, os, re, getopt, glob, numpy as np
from subprocess import call
from multiprocessing import Pool

## USER INPUT
PROCESSORS = 20
number_of_seqs = 5000
faa_dir = "/home/sam/FullCyc_metagenome/annotation/IMG/"
faa_file_prefix = "Ga0334612_proteins"

def signalP(fasta_chunk):
	print(fasta_chunk)
	os.system(' '.join([
		"signalp",
		"-fasta",
		fasta_chunk,
		"-org",
		"gram+",
		"-prefix",
		re.sub(".faa","gram_pos",fasta_chunk)
	]))

# Chunk Fasta
os.system(' '.join([
	"./fasta-splitter.pl",
	"--part-size",
	str(number_of_seqs),
	"--line-length 0 --measure count",
	faa_dir+faa_file_prefix+".faa"
	]))

# Grab File Names
file_list = []
for file in glob.glob(faa_file_prefix+".part*.faa"):
	file_list.append(file)
# Run Signal in Parallel
with Pool(processes=PROCESSORS) as pool:
	result = pool.map(signalP, file_list)
