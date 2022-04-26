#!/usr/bin/env python3


#python clean_up_input_fasta.py /Volumes/Seagate/Backup-MAC_HD2/proj_Alphonsine/Scripts_blast-positions/Blast-new-genome-pos/Genotype_F2-SNPs-raw_probe.fa /Volumes/Seagate/Backup-MAC_HD2/proj_Alphonsine/Scripts_blast-positions/Blast-new-genome-pos/Genotype_F2-SNPs-raw_probe-clean.fa
"""
python 5.0_clean_up_anyinput_fasta.py /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/PlantTFDBv5.0/Ath_pep.fas /home/gala0002/proj/proj_dir/REF/Arabidopsis-thaliana/PlantTFDBv5.0/Ath_pep_pars.fas
python 5.0_clean_up_anyinput_fasta.py /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/PlantTFDBv5.0/Gma_pep.fas /home/gala0002/proj/proj_dir/REF/Glycine-max_soybean/Glycine_Max_TF/PlantTFDBv5.0/Gma_pep_pars.fas
python 5.0_clean_up_anyinput_fasta.py /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/PlantTFDBv5.0/Mtr_pep.fas /home/gala0002/proj/proj_dir/REF/Medicago-truncatula/Medicago_Truncatula_TF/PlantTFDBv5.0/Mtr_pep_pars.fas


"""

import sys
import Bio

from Bio import SeqIO

input_fasta = sys.argv[1]
outputfl = open(sys.argv[2],"w")
for record in SeqIO.parse(input_fasta, "fasta"):
	#hdr=str(" [organism=Lepidium campestre]​")
	header_pars = record.description.split()
	header_pars2 = header_pars[0]+"|"+header_pars[2].split("|")[1]
	#print(header_pars2, sequence)
	fasta_entry = (">{}\n{}\n".format(header_pars2,record.seq))
	outputfl.write(str(fasta_entry))
outputfl.close()
