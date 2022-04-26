#!/bin/bash -l



########

#./2.0_de-novo-assembly_Trinity_all_sampleinfo.sh > 2.0_de-novo-assembly_Trinity_all_sampleinfo.sh.log 2>&1

#Trinity="/bioinfo/Trinity/trinityrnaseq-v2.11.0/Trinity"
TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"

INDIR="/home/gala0002/proj/proj_dir/1.0_sort-trim/"
OUTDIR="/home/gala0002/proj/proj_dir/2.0_trinity_out_all_sampleinfo"

nice -5 ${TRINITY_HOME}Trinity --seqType fq \
--samples_file /home/gala0002/proj/proj_dir/sampleinfo.txt --CPU 70 --max_memory 650G --output ${OUTDIR}

#----

${TRINITY_HOME}/util/TrinityStats.pl ${OUTDIR}/Trinity.fasta


#######################
#https://www.protocols.io/view/de-novo-transcriptome-assembly-workflow-ghebt3e?step=10
#####
# CD-HIT-EST
#####

WORKDIR="/home/gala0002/proj/proj_dir"
#for Trinity
cd $WORKDIR/
mkdir $WORKDIR/2.1_trinity_cdhit

cd-hit-est -i $WORKDIR/2.0_trinity_out_all_sampleinfo/Trinity.fasta -o $WORKDIR/2.1_trinity_cdhit/Trinity_clustered_0.95 -c 0.95 -n 8 -p 1 -g 1 -M 200000 -T 70 -d 40 1>$WORKDIR/2.1_trinity_cdhit/Trinity_cd-hit-est_0.95.log 2>$WORKDIR/2.1_trinity_cdhit/Trinity_cd-hit-est_0.95.err

cd-hit-est -i ../2.1.1_trinity_gene_splice_modeler/trinity_gene_splice_modeler.fasta -o $WORKDIR/2.1_trinity_cdhit/trinity_gene_splice_modeler_clustered_0.99 -c 0.99 -n 8 -p 1 -g 1 -M 200000 -T 70 -d 40 1>$WORKDIR/2.1_trinity_cdhit/trinity_gene_splice_modeler_cd-hit-est_0.99.log 2>$WORKDIR/2.1_trinity_cdhit/trinity_gene_splice_modeler_cd-hit-est_0.99.err

######

<<COMM1

#Download and install EvidentialGene tr2aacds
wget ftp://arthropods.eugenes.org/evigene.tar
tar xvf evigene.tar
cd evigene/
#add executables from scripts directory and subdirectories to PATH
echo export PATH=$PATH$( find $PROGRAMDIR/evigene/scripts/ -type d -printf 

######
#https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
#To install using ruby
cd $PROGRAMDIR
gem install transrate
transrate --install-deps all
#to install from source
cd $PROGRAMDIR
wget https://bintray.com/artifact/download/blahah/generic/transrate-1.0.2-linux-x86_64.tar.gz -O transrate.tar.gz
tar -zxvf transrate.tar.gz
cd transrate-1.0.2-linux-x86_64/bin
#add executables to a directory contained in PATH, or add current directory to PATH
echo export PATH=\$PATH:`pwd`\ >> ~/.bashrc && source ~/.bashrc
#Prior to run, edit /transrate/lib/app/lib/transrate/snap.rb to increase MCP
#edit the function “build_paired_cmd”, by finding the line “> cmd << 

COMM1

#/bioinfo/transrate-1.0.3-linux-x86_64

/bioinfo/transrate-1.0.3-linux-x86_64/transrate --assembly 2.0_trinity_out_all_sampleinfo/Trinity.fasta --threads 70 --output 2.1.2_trinity_transrate/transrate


########
# run BUSCO 
#######
#https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
# install
#conda create -n busco -c bioconda busco
#conda activate busco

#######

#Example:
#run_BUSCO.py -c 66 -o trinity_BUSCO_fabales_odb10 -i /home/gala0002/proj/proj_dir/2.0_trinity_out_all_sampleinfo/Trinity.fasta -l /bioinfo/BUSCO/busco/v1_lineage_sets/fabales_odb10 -m transcriptome

run_BUSCO.py -c 66 -o trinity_BUSCO_fabales_odb10 -i /home/gala0002/proj/proj_dir/2.0_trinity_out_all_sampleinfo/Trinity.fasta -l /bioinfo/BUSCO/busco/v1_lineage_sets/fabales_odb10 -m transcriptome

run_BUSCO.py -c 68 -o trinity_BUSCO_eudicots_odb10 -i /home/gala0002/proj/proj_dir/2.0_trinity_out_all_sampleinfo/Trinity.fasta -l /bioinfo/BUSCO/busco/v1_lineage_sets/eudicots_odb10 -m transcriptome


########
# run Trinity_gene_splice_modeler.py
#$TRINITY_HOME/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta Trinity.fasta
TRINITY_HOME="/bioinfo/trinityrnaseq-Trinity-v2.8.5/"


${TRINITY_HOME}/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
--trinity_fasta ./2.0_trinity_out_all_sampleinfo/Trinity.fasta --out_prefix 2.1.1_trinity_gene_splice_modeler



