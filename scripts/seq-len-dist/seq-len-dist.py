#!/usr/bin/python


# cat Trinity.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Trinity_seqlength-stats.txt

#conda activate py3.7

import Bio, pandas
from Bio import SeqIO

lengths = map(len, Bio.SeqIO.parse('/home/gala0002/proj/proj_dir/2.0_trinity_out_all_sampleinfo/Trinity.fasta', 'fasta'))


import matplotlib.pyplot as plt
pandas.Series(lengths).hist(color='gray', bins=100)
plt.savefig('/home/gala0002/proj/proj_dir/test_all/seq-len-hist.pdf')  # save the current figure
