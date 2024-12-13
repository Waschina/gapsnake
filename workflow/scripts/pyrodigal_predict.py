import os
import re
from Bio import SeqIO
import pyrodigal
from xopen import xopen
import gzip

genome = snakemake.input.genome
translation_table = int(snakemake.params.translation_table)

# read contigs
sequences = []
with xopen(genome, mode='r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sequences.append(str(record.seq))

# train pyrodigal
orf_finder = pyrodigal.GeneFinder(closed=False)
orf_finder.train(*sequences, translation_table=translation_table)

# predict genes
f_faa = open(os.path.splitext(re.sub('\\.gz$','',genome))[0]+'.faa', 'w+')

with xopen(genome, mode='r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        genes = orf_finder.find_genes(bytes(record.seq))
        with open(f_faa.name, "a") as dst:
            genes.write_translations(dst, sequence_id=record.id, width=80)

# gzip protein fasta
with open(f_faa.name, 'rb') as f_in, gzip.open(snakemake.output.prot, 'wb') as f_out:
    f_out.writelines(f_in)
os.remove(f_faa.name)

