import re
import pandas as pd
import plotly.graph_objects as go
import subprocess

# this is a bit ugly, should find a way to do it more efficiently
R1 = len([1 for line in open(snakemake.input.R1) if re.match('^\+$', line)])
R2 = len([1 for line in open(snakemake.input.R2) if re.match('^\+$', line)])
paired = len([1 for line in open(snakemake.input.R1R2) if re.match('^\+$', line)])
good_ali = len([1 for line in open(snakemake.input.good) if re.match('^\+$', line)])
good_ali_p = round(float(good_ali)/float(paired)*100,1)
bad_ali = len([1 for line in open(snakemake.input.bad) if re.match('^\+$', line)])
bad_ali_p = round(float(bad_ali)/float(paired)*100,1)
demult = len([1 for line in open(snakemake.input.demult) if re.match('^>', line)])
demult_p = round(float(demult)/float(good_ali)*100,1)
unassigned = len([1 for line in open(snakemake.input.unass) if re.match('^>', line)])
unassigned_p = round(float(unassigned)/float(good_ali)*100,1)

df = pd.DataFrame(columns = ['file', 'reads', 'perc_kept_lost'],
	index = ['R1', 'R2', 'paired', 'good_ali', 'bad_ali', 'demultiplexed', 'unassigned'])

df['file'] = ['R1', 'R2', 'paired', 'good_ali', 'bad_ali', 'demultiplexed', 'unassigned']
df['reads'] = [R1, R2, paired, good_ali, bad_ali, demult, unassigned]
df['kept_lost'] = ['-', '-', '-', good_ali_p, bad_ali_p, demult_p, unassigned_p]

fig = go.Figure(data = [go.Table(
	header = dict(values = ['', '# reads/seq', '% kept/lost'],
		line_color='black',
		fill_color = 'grey',
		font=dict(color='white')
		),
	cells = dict(values = [df.file, df.reads, df.kept_lost],
		line_color='black',
		fill_color = 'white'))
])

fig.write_image(snakemake.output.svg1)

df.to_csv(snakemake.output.tab1, sep='\t')


def count_all_seq(obifasta):
	n = 0
	file = open(obifasta, "r")
	for line in file:
		if re.search("^>", line):
			a = line.split("; ")
			for value in a:
				if 'count' in value:
					b = value.split("=")
					n += int(b[1])
	return(n)


derepl_uniq = len([1 for line in open(snakemake.input.derepl) if re.match('^>', line)])
derepl_seq = count_all_seq(snakemake.input.derepl)

basicfilt_uniq = len([1 for line in open(snakemake.input.filt) if re.match('^>', line)])
basicfilt_seq = count_all_seq(snakemake.input.filt)
basicfilt_uniq_p = round(float(basicfilt_uniq)/float(derepl_uniq)*100,1)
basicfilt_seq_p = round(float(basicfilt_seq)/float(derepl_seq)*100,1)

clust_uniq = len([1 for line in open(snakemake.input.clust) if re.match('^>', line)])
clust_seq = count_all_seq(snakemake.input.clust)

df = pd.DataFrame(columns = ['file', 'uniq_seq', 'perc_kept_uniq', 'seq', 'perc_kept_seq', 'motus'],
	index = ['dereplicated', 'basicfilt', 'clustering'])

df['file'] = ['dereplicated', 'basicfilt', 'clustering']
df['uniq_seq'] = [derepl_uniq, basicfilt_uniq, '-']
df['kept_uniq'] = ['-', '-', basicfilt_uniq_p]
df['seq'] = [derepl_seq, basicfilt_seq, clust_seq]
df['kept_seq'] = ['-', basicfilt_seq_p, '-']
df['motus'] = ['-', '-', clust_uniq]

fig = go.Figure(data = [go.Table(
	header = dict(values = ['', '# uniq seq', '% kept uniq', '# seq', '% kept seq', '# MOTUs'],
		line_color='black',
		fill_color = 'grey',
		font=dict(color='white')
		),
	cells = dict(values = [df.file, df.uniq_seq, df.kept_uniq, df.seq, df.kept_seq, df.motus],
		line_color='black',
		fill_color = 'white'))
])

fig.write_image(snakemake.output.svg2)

df.to_csv(snakemake.output.tab2, sep='\t')

