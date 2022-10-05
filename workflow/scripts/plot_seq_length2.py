import seaborn as sns
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


sizes = [len(rec) for rec in SeqIO.parse(snakemake.input.cs, 'fasta')]
# plt.hist(sizes, color = 'green')
# plt.title('Histogram of sequences length')
# plt.xlabel('Sequence length (bp)')
# plt.ylabel('Count')
# plt.savefig(snakemake.output.csl)
# plt.clf()

sns.set_style("whitegrid")
p = sns.displot(sizes)
p.set(title = 'Histogram of sequences length after clustering', xlabel = "Sequences length")
plt.savefig(snakemake.output.csl)
plt.clf()