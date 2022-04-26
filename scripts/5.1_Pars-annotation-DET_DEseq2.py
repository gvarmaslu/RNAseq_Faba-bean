#!/usr/bin/python


"""
#Script to Pars vcf...

#####------Inputs-------
# Run-RNA-workflow-7.0_Pars-annotation-DET.py INPUT=Inputfile OUTPUT=Outputfile DIR=Directory-fullpath

#python 5.1_Pars-annotation-DET_DEseq2.py F-EMI_vs_F-EMII.DESeq2.DE_results.P5e-2_C1.DE.subset_head F-EMI_vs_F-EMII.DESeq2.DE_results.P5e-2_C1.DE.subset_head_anno.tsv /home/gala0002/proj/proj_dir/DESeq2_genes_F-EMI_vs_F-EMII/


#############################################------
#cut -d";" --output-delimiter=$'\t' -f3,4,5,7,8,11 | 

truncate -s 0 Thaliana_PIDs_anno

####------


"""

import sys
import re
import subprocess
from subprocess import *
from subprocess import call

"""
#Query_ID, Subject_ID, %Identity, Alignment_length, Mismatches, Gap_opens, Q.start, Q.end, S.start, S.end, E-value, Bit-score
#Query_ID	Subject_ID	%Identity	Alignment_length	Mismatches	Gap_opens	Q.start	Q.end	S.start	S.end	E-value	Bit-score
"""


class fileHandler:
	def __init__(self):
		self.data = []
		#print "Calling fileHandler constructor"
	def open_file(self,readfl):
		self.rfile = open(readfl,'r').readlines()
		return self.rfile
	def write_file(self,writefl):
		self.wfile = open(writefl,'w')
		return self.wfile

class SearchDB(fileHandler):

	def __init__(self):
		self.data = []
		from collections import defaultdict
		self.ident_ranges_HMBM = defaultdict(list)

	def Search_DBs(self,readfl1,outfl1,workdir):
		"""
		Calling Search local DB
		"""
		def fun_unquote(string):
			''' Decode escaped characters in URL: Returns the string without escaped/ASCII characters'''
			import urllib.parse
			return urllib.parse.unquote(string)

		def blastout_pars(AnnoPars):
			GeneID2 = AnnoPars[-2].strip().split()[1].replace(':', ' ').replace('=', ' ').replace(',', ' ').split()
			if len(GeneID2) > 2:
				GeneID2.sort()
				GeneID2 = GeneID2[0]+"@"+GeneID2[1]
			else:
				GeneID2 = GeneID2[1]
			product= AnnoPars[-1].split("=")[1]
			qstart_end=str("Q:"+"-".join(AnnoPars[6:8]))
			sstart_end=str("S:"+"-".join(AnnoPars[8:10]))
			AnnoPars_merge= str(qstart_end+","+sstart_end+";"+AnnoPars[2]+"%ID;E:"+AnnoPars[10])	
			grepout = str(AnnoPars[1]+"\t"+AnnoPars_merge+"\t"+GeneID2+"\t"+fun_unquote(product))
			return grepout

		def srchdb_blastout(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\t")
				grepout = blastout_pars(AnnoPars)
			except:
				False
				grepout = str("."+"\t"+"."+"\t"+"."+"\t"+".")
			return grepout

		def srchdb_blastout_transdecoder(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep '"+str(GName)+"' "+str(FPATH)+" | sort -k9,10"
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\n")
				if len(AnnoPars) > 1:
					pepIDs=[];pepDisc=[]
					for pp in AnnoPars:
						ppep = pp.split("\t")
						BP_ppep = blastout_pars(ppep)
						BP_ppep = BP_ppep.replace("\t",";")
						pepIDs.append(ppep[0])
						pepDisc.append(BP_ppep)
					pepIDs_jn="@".join(pepIDs)
					pepDisc_jn="@".join(pepDisc)
					grepout = str(pepIDs_jn+"\t"+pepDisc_jn)
				elif (AnnoPars != ['']):
					ppep = AnnoPars[0].split("\t")
					BP_ppep= blastout_pars(ppep)
					pepIDs_jn= ppep[0]
					pepDisc_jn= BP_ppep.replace("\t",";")
					grepout = str(pepIDs_jn+"\t"+pepDisc_jn)
				else:
					grepout = str("."+"\t"+".")
			except:
				False
				grepout = str("."+"\t"+".")
			return grepout

		def srchdb_blastout_TFDB(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 '"+str(GName)+"' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\t")
				qstart_end=str("Q:"+"-".join(AnnoPars[6:8]))
				sstart_end=str("S:"+"-".join(AnnoPars[8:10]))
				AnnoPars_merge= str(qstart_end+","+sstart_end+";"+AnnoPars[2]+"%ID;E:"+AnnoPars[10])
				grepout = str(AnnoPars[1]+"\t"+AnnoPars_merge)
			except:
				False
				grepout = str("."+"\t"+".")
			return grepout

		def srchdb_Trinotate(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep -w '"+str(GName)+"' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\n")
				if len(AnnoPars) > 1:
					geneDisc=[];pepIDs=[];pepDisc=[]
					for pp in AnnoPars:
						ppep = pp.split("\t")
						pgene = ppep[2]
						ppep_disc = ";".join(ppep[4:])
						geneDisc.append(pgene)
						pepIDs.append(ppep[4])
						pepDisc.append(ppep_disc)
					pepIDs_jn="@".join(pepIDs)
					pepDisc_jn="@".join(pepDisc)
					grepout = str(geneDisc[0]+"\t"+pepIDs_jn+"\t"+pepDisc_jn)
				elif len(AnnoPars) == 1:
					ppep = AnnoPars[0].split("\t")
					ppep_disc = ";".join(ppep[4:])
					grepout = str(ppep[2]+"\t"+ppep[4]+"\t"+ppep_disc)
				else:
					grepout = str("."+"\t"+"."+"\t"+".")
			except:
				False
				grepout = str("."+"\t"+"."+"\t"+".")
			return grepout

		def srchdb_Trinotate_GO(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep -m 1 -w '"+str(GName)+"' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\t")
				grepout = str(AnnoPars[1])
			except:
				False
				grepout = str(".")
			return grepout

		def srchdb_blastout_all_transdecoder(GName,FPATH):
			try:
				True
				cmdFls1 = "LANG=C grep '"+str(GName)+"' "+str(FPATH)+""
				cmdFls2 =  subprocess.check_output(cmdFls1, shell=True)
				AnnoPars = cmdFls2.strip().decode().split("\n")
				if len(AnnoPars) > 1:
					pepIDs=[];pepDisc=[]
					for pp in AnnoPars:
						ppep = pp.split("\t")
						BP_ppep = blastout_pars(ppep)
						BP_ppep = BP_ppep.replace("\t",";")
						pepIDs.append(ppep[0])
						pepDisc.append(BP_ppep)
					pepIDs_jn=pepIDs[0]
					pepDisc_jn=pepDisc[0]
					grepout = str(pepIDs_jn+"\t"+pepDisc_jn)
				elif (AnnoPars != ['']):
					ppep = AnnoPars[0].split("\t")
					BP_ppep= blastout_pars(ppep)
					pepIDs_jn= ppep[0]
					pepDisc_jn= BP_ppep.replace("\t",";")
					grepout = str(pepIDs_jn+"\t"+pepDisc_jn)
				else:
					grepout = str("."+"\t"+".")
			except:
				False
				grepout = str("."+"\t"+".")
			return grepout
		with open(workdir+readfl1,'r') as f1, open(workdir+outfl1,'w') as output:
			first_line0 = f1.readline().strip().split("\t")
			first_lines = f1.readlines()
			#output.write(str("Quinoa_mRNA-ID"+"\t"+str("\t".join(first_line0))+"\t"+"Quinoa-ASM168347v1-CDS-ID	Protein-ID	Gene_ID	Gene_Name	Gene_Pos	Gene-Ann1	Gene-Ann2	Gene-Ann3	Gene-Ann4	Gene-Ann5	Target_AT_PID_HMMER	Target_AT_E-value_HMMER	Target_AT_Score_HMMER	Target_AT_Description_HMMER	Target_AT_PID_NCBIBLAST	Target_AT_E-value_NCBIBLAST	Target_AT_Score_NCBIBLAST	Target_AT_GeneID_NCBIBLAST	Target_AT_Description_NCBIBLAST	PlantTFDBv5.0_Ext-Cqu-cds_TFfamID	PlantTFDBv5.0_Evalue	PlantTFDBv5.0_Identity_Score	ITAKv18.12_Ext-Cqu-cds_TFfamID	ITAKv18.12_Evalue	ITAKv18.12_Identity_Score	Target_Cq.PI614886_PID_NCBIBLAST	Target_Cq.PI614886_E-value_NCBIBLAST	Target_Cq.PI614886_Score_NCBIBLAST"+"\n"))
			#### Insert headers
			Geneanno_blasthead= str("Q=Query:start-end,S=Subject:start-end;%ID=%Identity;E=E-value	GeneID	product")
			Geneanno_outhead= str("Subject:Medicago truncatula	"+Geneanno_blasthead+"\t"+"Subject:Pisum sativum (Pea)	"+Geneanno_blasthead+"\t"+"Subject:Glycine max (Soybean)	"+Geneanno_blasthead+"\t"+"Subject:Arabidopsis thaliana	"+Geneanno_blasthead)
			####
			TransD_blasthead= str("IsoformID1;Q=Query:start-end,S=Subject:start-end;%ID=%Identity;E=E-value;GeneID;product @ IsoformID2;Q=Query:start-end,S=Subject:start-end;%ID=%Identity;E=E-value;GeneID;product")
			TransD_outhead= str("Medicago truncatula (IsoformID1@IsoformID2)	"+TransD_blasthead+"\t"+"Pisum sativum (Pea) (IsoformID1@IsoformID2)	"+TransD_blasthead+"\t"+"Glycine max (Soybean) (IsoformID1@IsoformID2)	"+TransD_blasthead+"\t"+"Arabidopsis thaliana (IsoformID1@IsoformID2)	"+TransD_blasthead)
			####
			TF_blasthead = str("Q=Query:start-end,S=Subject:start-end;%ID=%Identity;E=E-value")
			TF_outhead= str("Medicago truncatula (PlantTFDB)	"+TF_blasthead+"\t"+"Medicago truncatula (iTAK_MT)	"+TF_blasthead+"\t"+"Glycine max (Soybean) (PlantTFDB)	"+TF_blasthead+"\t"+"Glycine max (Soybean) (iTAK_MT)	"+TF_blasthead+"\t"+"Arabidopsis thaliana (PlantTFDB)	"+TF_blasthead+"\t"+"Arabidopsis thaliana (iTAK_MT)	"+TF_blasthead)
			####
			Trino_blasthead= str("@IsoformID1;prot_id;prot_coords;sprot_Top_BLASTP_hit;Pfam;SignalP;TmHMM;eggnog;Kegg;gene_ontology_BLASTX;gene_ontology_BLASTP;gene_ontology_Pfam;transcript;peptide@")
			Trino_outhead= str("sprot_Top_BLASTX_hit(Gene_Description Swiss_Prot)	sprot_Top_BLASTP_hit(IsoformID1@IsoformID2)	"+Trino_blasthead+"\t"+"Gene Ontology IDs")
			####
			"""
			TransD2_blasthead= str("IsoformID1;Q=Query:start-end,S=Subject:start-end;%ID=%Identity;E=E-value;GeneID;product")
			TransD2_outhead= str("Medicago truncatula (IsoformID1)	"+TransD2_blasthead+"\t"+"Pisum sativum (Pea) (IsoformID1)	"+TransD2_blasthead+"\t"+"Glycine max (Soybean) (IsoformID1)	"+TransD2_blasthead+"\t"+"Arabidopsis thaliana (IsoformID1)	"+TransD2_blasthead)
			#+"\t"+TransD2_outhead+
			"""
			####
			output.write(str("GeneID_TRINITY"+"\t"+str("\t".join(first_line0))+"\t"+Geneanno_outhead+"\t"+TransD_outhead+"\t"+TF_outhead+"\t"+Trino_outhead+"\n"))
			####
			for lns in first_lines:
				lns_sp =  lns.strip().split("\t")
				lns_sp2 =  lns_sp[0]
				#print(lns_sp2)
				DIR="/home/gala0002/proj/proj_dir/4.2_blastout/"
				FPATH1=DIR+"trinity_gene_splice_modeler.fasta_Medicago-truncatula.blastx.outfmt6_anno.tsv"
				FPATH2=DIR+"trinity_gene_splice_modeler.fasta_Pisum-sativum.blastx.outfmt6_anno.tsv"
				FPATH3=DIR+"trinity_gene_splice_modeler.fasta_Glycine-max.blastx.outfmt6_anno.tsv"
				FPATH4=DIR+"trinity_gene_splice_modeler.fasta_Arabidopsis-thaliana.blastx.outfmt6_anno.tsv"

				#####
				# Gene splice modeler Gene annotation 
				#####
				#Medicago truncatula
				Anno_out1 = srchdb_blastout(lns_sp2,FPATH1)
				#Pisum sativum (Pea)
				Anno_out2 = srchdb_blastout(lns_sp2,FPATH2)
				#Glycine max (Soybean)
				Anno_out3 = srchdb_blastout(lns_sp2,FPATH3)
				#Arabidopsis thaliana
				Anno_out4 = srchdb_blastout(lns_sp2,FPATH4)
				#print(Anno_out1+"\t"+Anno_out2+"\t"+Anno_out3+"\t"+Anno_out4)
				#print(Anno_out1)
				Geneanno_out= str(Anno_out1+"\t"+Anno_out2+"\t"+Anno_out3+"\t"+Anno_out4)

				#####
				# Transdecoder Isoforms annotation 
				#####
				DIR="/home/gala0002/proj/proj_dir/4.2_blastout/"
				FPATH_TD1=DIR+"trinity_gene_splice_modeler.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6_anno.tsv"
				FPATH_TD2=DIR+"trinity_gene_splice_modeler.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6_anno.tsv"
				FPATH_TD3=DIR+"trinity_gene_splice_modeler.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6_anno.tsv"
				FPATH_TD4=DIR+"trinity_gene_splice_modeler.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6_anno.tsv"
				#Medicago truncatula
				Anno_out_TD1 = srchdb_blastout_transdecoder(lns_sp2,FPATH_TD1)
				#Pisum sativum (Pea)
				Anno_out_TD2 = srchdb_blastout_transdecoder(lns_sp2,FPATH_TD2)
				#Glycine max (Soybean)
				Anno_out_TD3 = srchdb_blastout_transdecoder(lns_sp2,FPATH_TD3)
				#Arabidopsis thaliana
				Anno_out_TD4 = srchdb_blastout_transdecoder(lns_sp2,FPATH_TD4)
				#print(Anno_out_TD1+"\t"+Anno_out_TD2+"\t"+Anno_out_TD3+"\t"+Anno_out_TD4)
				#print(Anno_out_TD1)
				TransD_out= str(Anno_out_TD1+"\t"+Anno_out_TD2+"\t"+Anno_out_TD3+"\t"+Anno_out_TD4)


				#####
				# Transcription Factor (TF) database annotation 
				#####
				#PlantTFDBv5.0
				DIR="/home/gala0002/proj/proj_dir/4.2_blastout/"
				FPATH_PlantTFDB_MT=DIR+"trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Medicago-truncatula.blastx.outfmt6"
				FPATH_PlantTFDB_GM=DIR+"trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Glycine-max.blastx.outfmt6"
				FPATH_PlantTFDB_AT=DIR+"trinity_gene_splice_modeler.fasta_PlantTFDBv5.0_Arabidopsis-thaliana.blastx.outfmt6"
				#iTAKv18.12
				DIR="/home/gala0002/proj/proj_dir/4.2_blastout/"
				FPATH_iTAK_MT=DIR+"trinity_gene_splice_modeler.fasta_iTAKv18.12_Medicago-truncatula.blastx.outfmt6"
				FPATH_iTAK_GM=DIR+"trinity_gene_splice_modeler.fasta_iTAKv18.12_Glycine-max.blastx.outfmt6"
				FPATH_iTAK_AT=DIR+"trinity_gene_splice_modeler.fasta_iTAKv18.12_Arabidopsis-thaliana.blastx.outfmt6"
				#####
				#Medicago truncatula
				Anno_out_TF_MT1 = srchdb_blastout_TFDB(lns_sp2,FPATH_PlantTFDB_MT)
				Anno_out_TF_MT2 = srchdb_blastout_TFDB(lns_sp2,FPATH_iTAK_MT)
				#Glycine max (Soybean)
				Anno_out_TF_GM1 = srchdb_blastout_TFDB(lns_sp2,FPATH_PlantTFDB_GM)
				Anno_out_TF_GM2 = srchdb_blastout_TFDB(lns_sp2,FPATH_iTAK_GM)
				#Arabidopsis thaliana
				Anno_out_TF_AT1 = srchdb_blastout_TFDB(lns_sp2,FPATH_PlantTFDB_AT)
				Anno_out_TF_AT2 = srchdb_blastout_TFDB(lns_sp2,FPATH_iTAK_AT)
				#print(Anno_out_TF_MT1+"\t"+Anno_out_TF_MT2+"\t"+Anno_out_TF_GM1+"\t"+Anno_out_TF_GM2+"\t"+Anno_out_TF_AT1+"\t"+Anno_out_TF_AT2)
				TF_out= str(Anno_out_TF_MT1+"\t"+Anno_out_TF_MT2+"\t"+Anno_out_TF_GM1+"\t"+Anno_out_TF_GM2+"\t"+Anno_out_TF_AT1+"\t"+Anno_out_TF_AT2)

				#####
				# Trinotate, GO annotation 
				#####
				DIR="/home/gala0002/proj/proj_dir/4.1.2_Trinotate/"
				FPATH_Trino=DIR+"trinotate_annotation_report.xls"
				FPATH_Trino_GO=DIR+"trinotate_annotation_report_GO.txt"
				#####
				Anno_out_Trino = srchdb_Trinotate(lns_sp2,FPATH_Trino)
				Anno_out_Trino_GO = srchdb_Trinotate_GO(lns_sp2,FPATH_Trino_GO)
				#print(Anno_out_Trino+"\t"+Anno_out_Trino_GO)
				Trino_out=str(Anno_out_Trino+"\t"+Anno_out_Trino_GO)

				#####
				# Transdecoder Isoforms annotation (All Trinity transcripts)
				#####
				'''
				DIR="/home/gala0002/proj/proj_dir/4.2_blastout/"
				FPATH2_TD1=DIR+"Trinity.fasta.transdecoder.pep_Medicago-truncatula.blastp.outfmt6_anno.tsv"
				FPATH2_TD2=DIR+"Trinity.fasta.transdecoder.pep_Pisum-sativum.blastp.outfmt6_anno.tsv"
				FPATH2_TD3=DIR+"Trinity.fasta.transdecoder.pep_Glycine-max.blastp.outfmt6_anno.tsv"
				FPATH2_TD4=DIR+"Trinity.fasta.transdecoder.pep_Arabidopsis-thaliana.blastp.outfmt6_anno.tsv"
				#Medicago truncatula
				Anno_out2_TD1 = srchdb_blastout_all_transdecoder(lns_sp2,FPATH2_TD1)
				#Pisum sativum (Pea)
				Anno_out2_TD2 = srchdb_blastout_all_transdecoder(lns_sp2,FPATH2_TD2)
				#Glycine max (Soybean)
				Anno_out2_TD3 = srchdb_blastout_all_transdecoder(lns_sp2,FPATH2_TD3)
				#Arabidopsis thaliana
				Anno_out2_TD4 = srchdb_blastout_all_transdecoder(lns_sp2,FPATH2_TD4)
				#print(Anno_out_TD1+"\t"+Anno_out_TD2+"\t"+Anno_out_TD3+"\t"+Anno_out_TD4)
				#print(Anno_out2_TD1)
				TransD2_out= str(Anno_out2_TD1+"\t"+Anno_out2_TD2+"\t"+Anno_out2_TD3+"\t"+Anno_out2_TD4)
				Anno_outall = str("\t".join(lns_sp)+"\t"+Geneanno_out+"\t"+TransD_out+"\t"+TF_out+"\t"+Trino_out+"\t"+TransD2_out+"\n")
				'''
				###########
				Anno_outall = str("\t".join(lns_sp)+"\t"+Geneanno_out+"\t"+TransD_out+"\t"+TF_out+"\t"+Trino_out+"\n")
				output.write(str(Anno_outall))
				
		print("Done search...")
		return None

clF1 = SearchDB().Search_DBs(sys.argv[1],sys.argv[2],sys.argv[3])





