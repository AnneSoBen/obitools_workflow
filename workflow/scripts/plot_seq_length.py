import seaborn as sns
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

length_paired = [len(rec) for rec in SeqIO.parse(snakemake.input.ps, 'fastq-illumina')]
length_alifilt = [len(rec) for rec in SeqIO.parse(snakemake.input.afs, 'fastq-illumina')]
length_demult = [len(rec) for rec in SeqIO.parse(snakemake.input.ds, 'fasta')]
s = length_paired + length_alifilt + length_demult
x = np.array(["paired", "alifilt", "demultiplexed"])
y = (len(length_paired), len(length_alifilt), len(length_demult))
r = np.repeat(x, y, axis=0)

df = pd.DataFrame(columns = ['seq_length', 'step'])
df['seq_length'] = s
df['step'] = r

sns.set_style("whitegrid")
p = sns.displot(df, x = "seq_length", hue = "step")
p.set(title = 'Histogram of sequences length', xlabel = "Sequences length")
plt.savefig(snakemake.output.sl)
plt.clf()
