#!/bin/bash -l



########

#./2.0_de-novo-assembly_Trinity_all.sh > 2.0_de-novo-assembly_Trinity_all.sh.log 2>&1

#Trinity="/bioinfo/Trinity/trinityrnaseq-v2.11.0/Trinity"
Trinity="/bioinfo/trinityrnaseq-Trinity-v2.8.5/Trinity"

INDIR="/home/gala0002/proj/proj_dir/1.0_sort-trim/"


nice -5 ${Trinity} --seqType fq \
--left ${INDIR}F1_EMIII_lib435236_7094_4/F1_EMIII_lib435236_7094_4-sort-trim_1.fq.gz,${INDIR}F1_EMII_lib435231_7093_2/F1_EMII_lib435231_7093_2-sort-trim_1.fq.gz,${INDIR}F1_EMI_lib435225_7093_2/F1_EMI_lib435225_7093_2-sort-trim_1.fq.gz,${INDIR}F1_ESII_lib435246_7093_3/F1_ESII_lib435246_7093_3-sort-trim_1.fq.gz,${INDIR}F1_ESI_lib435241_7093_2/F1_ESI_lib435241_7093_2-sort-trim_1.fq.gz,${INDIR}F1_PS_lib435254_7093_3/F1_PS_lib435254_7093_3-sort-trim_1.fq.gz,${INDIR}F2_EMI_lib435226_7093_2/F2_EMI_lib435226_7093_2-sort-trim_1.fq.gz,${INDIR}F2_ESII_lib435247_7093_3/F2_ESII_lib435247_7093_3-sort-trim_1.fq.gz,${INDIR}F2_ESI_lib435242_7093_2/F2_ESI_lib435242_7093_2-sort-trim_1.fq.gz,${INDIR}F2_PS_lib435255_7093_3/F2_PS_lib435255_7093_3-sort-trim_1.fq.gz,${INDIR}F3_EMIII_lib435237_7093_2/F3_EMIII_lib435237_7093_2-sort-trim_1.fq.gz,${INDIR}F3_ESII_lib435248_7093_3/F3_ESII_lib435248_7093_3-sort-trim_1.fq.gz,${INDIR}F4_PS_lib435256_7093_3/F4_PS_lib435256_7093_3-sort-trim_1.fq.gz,${INDIR}F5_EMIII_lib435238_7093_2/F5_EMIII_lib435238_7093_2-sort-trim_1.fq.gz,${INDIR}F5_EMII_lib435232_7093_2/F5_EMII_lib435232_7093_2-sort-trim_1.fq.gz,${INDIR}F5_PII_lib435252_7093_3/F5_PII_lib435252_7093_3-sort-trim_1.fq.gz,${INDIR}F6_EMI_lib435227_7093_2/F6_EMI_lib435227_7093_2-sort-trim_1.fq.gz,${INDIR}F6_ESI_lib435243_7093_3/F6_ESI_lib435243_7093_3-sort-trim_1.fq.gz,${INDIR}T3_EMI_lib435228_7093_2/T3_EMI_lib435228_7093_2-sort-trim_1.fq.gz,${INDIR}T4_EMII_lib441948_7094_4/T4_EMII_lib441948_7094_4-sort-trim_1.fq.gz,${INDIR}T4_EMI_lib435229_7093_2/T4_EMI_lib435229_7093_2-sort-trim_1.fq.gz,${INDIR}T4_ESII_lib435249_7090_3/T4_ESII_lib435249_7090_3-sort-trim_1.fq.gz,${INDIR}T4_PS_lib435257_7093_3/T4_PS_lib435257_7093_3-sort-trim_1.fq.gz,${INDIR}T5_EMIII_lib435239_7093_2/T5_EMIII_lib435239_7093_2-sort-trim_1.fq.gz,${INDIR}T5_EMII_lib435234_7093_2/T5_EMII_lib435234_7093_2-sort-trim_1.fq.gz,${INDIR}T5_EMI_lib435230_7093_2/T5_EMI_lib435230_7093_2-sort-trim_1.fq.gz,${INDIR}T5_ESII_lib435250_7093_3/T5_ESII_lib435250_7093_3-sort-trim_1.fq.gz,${INDIR}T5_ESI_lib435244_7093_3/T5_ESI_lib435244_7093_3-sort-trim_1.fq.gz,${INDIR}T5_PII_lib435253_7093_3/T5_PII_lib435253_7093_3-sort-trim_1.fq.gz,${INDIR}T6_EMII_lib435235_7093_2/T6_EMII_lib435235_7093_2-sort-trim_1.fq.gz,${INDIR}T6_ESII_lib435251_7093_3/T6_ESII_lib435251_7093_3-sort-trim_1.fq.gz,${INDIR}T6_ESI_lib435245_7093_3/T6_ESI_lib435245_7093_3-sort-trim_1.fq.gz,${INDIR}T7_PS_lib435258_7093_3/T7_PS_lib435258_7093_3-sort-trim_1.fq.gz,${INDIR}T8_EMIII_lib435240_7093_2/T8_EMIII_lib435240_7093_2-sort-trim_1.fq.gz,${INDIR}T8_PS_lib435259_7093_3/T8_PS_lib435259_7093_3-sort-trim_1.fq.gz \
--right ${INDIR}F1_EMIII_lib435236_7094_4/F1_EMIII_lib435236_7094_4-sort-trim_2.fq.gz,${INDIR}F1_EMII_lib435231_7093_2/F1_EMII_lib435231_7093_2-sort-trim_2.fq.gz,${INDIR}F1_EMI_lib435225_7093_2/F1_EMI_lib435225_7093_2-sort-trim_2.fq.gz,${INDIR}F1_ESII_lib435246_7093_3/F1_ESII_lib435246_7093_3-sort-trim_2.fq.gz,${INDIR}F1_ESI_lib435241_7093_2/F1_ESI_lib435241_7093_2-sort-trim_2.fq.gz,${INDIR}F1_PS_lib435254_7093_3/F1_PS_lib435254_7093_3-sort-trim_2.fq.gz,${INDIR}F2_EMI_lib435226_7093_2/F2_EMI_lib435226_7093_2-sort-trim_2.fq.gz,${INDIR}F2_ESII_lib435247_7093_3/F2_ESII_lib435247_7093_3-sort-trim_2.fq.gz,${INDIR}F2_ESI_lib435242_7093_2/F2_ESI_lib435242_7093_2-sort-trim_2.fq.gz,${INDIR}F2_PS_lib435255_7093_3/F2_PS_lib435255_7093_3-sort-trim_2.fq.gz,${INDIR}F3_EMIII_lib435237_7093_2/F3_EMIII_lib435237_7093_2-sort-trim_2.fq.gz,${INDIR}F3_ESII_lib435248_7093_3/F3_ESII_lib435248_7093_3-sort-trim_2.fq.gz,${INDIR}F4_PS_lib435256_7093_3/F4_PS_lib435256_7093_3-sort-trim_2.fq.gz,${INDIR}F5_EMIII_lib435238_7093_2/F5_EMIII_lib435238_7093_2-sort-trim_2.fq.gz,${INDIR}F5_EMII_lib435232_7093_2/F5_EMII_lib435232_7093_2-sort-trim_2.fq.gz,${INDIR}F5_PII_lib435252_7093_3/F5_PII_lib435252_7093_3-sort-trim_2.fq.gz,${INDIR}F6_EMI_lib435227_7093_2/F6_EMI_lib435227_7093_2-sort-trim_2.fq.gz,${INDIR}F6_ESI_lib435243_7093_3/F6_ESI_lib435243_7093_3-sort-trim_2.fq.gz,${INDIR}T3_EMI_lib435228_7093_2/T3_EMI_lib435228_7093_2-sort-trim_2.fq.gz,${INDIR}T4_EMII_lib441948_7094_4/T4_EMII_lib441948_7094_4-sort-trim_2.fq.gz,${INDIR}T4_EMI_lib435229_7093_2/T4_EMI_lib435229_7093_2-sort-trim_2.fq.gz,${INDIR}T4_ESII_lib435249_7090_3/T4_ESII_lib435249_7090_3-sort-trim_2.fq.gz,${INDIR}T4_PS_lib435257_7093_3/T4_PS_lib435257_7093_3-sort-trim_2.fq.gz,${INDIR}T5_EMIII_lib435239_7093_2/T5_EMIII_lib435239_7093_2-sort-trim_2.fq.gz,${INDIR}T5_EMII_lib435234_7093_2/T5_EMII_lib435234_7093_2-sort-trim_2.fq.gz,${INDIR}T5_EMI_lib435230_7093_2/T5_EMI_lib435230_7093_2-sort-trim_2.fq.gz,${INDIR}T5_ESII_lib435250_7093_3/T5_ESII_lib435250_7093_3-sort-trim_2.fq.gz,${INDIR}T5_ESI_lib435244_7093_3/T5_ESI_lib435244_7093_3-sort-trim_2.fq.gz,${INDIR}T5_PII_lib435253_7093_3/T5_PII_lib435253_7093_3-sort-trim_2.fq.gz,${INDIR}T6_EMII_lib435235_7093_2/T6_EMII_lib435235_7093_2-sort-trim_2.fq.gz,${INDIR}T6_ESII_lib435251_7093_3/T6_ESII_lib435251_7093_3-sort-trim_2.fq.gz,${INDIR}T6_ESI_lib435245_7093_3/T6_ESI_lib435245_7093_3-sort-trim_2.fq.gz,${INDIR}T7_PS_lib435258_7093_3/T7_PS_lib435258_7093_3-sort-trim_2.fq.gz,${INDIR}T8_EMIII_lib435240_7093_2/T8_EMIII_lib435240_7093_2-sort-trim_2.fq.gz,${INDIR}T8_PS_lib435259_7093_3/T8_PS_lib435259_7093_3-sort-trim_2.fq.gz \
--CPU 70 --max_memory 650G --output /home/gala0002/proj/proj_dir/2.0_trinity_out_all


####


<<COMM

##### JUST F1

nice -5 ${Trinity} --seqType fq \
--left ${INDIR}F1_EMIII_lib435236_7094_4/F1_EMIII_lib435236_7094_4-sort-trim_1.fq.gz,${INDIR}F1_EMII_lib435231_7093_2/F1_EMII_lib435231_7093_2-sort-trim_1.fq.gz,${INDIR}F1_EMI_lib435225_7093_2/F1_EMI_lib435225_7093_2-sort-trim_1.fq.gz,${INDIR}F1_ESII_lib435246_7093_3/F1_ESII_lib435246_7093_3-sort-trim_1.fq.gz,${INDIR}F1_ESI_lib435241_7093_2/F1_ESI_lib435241_7093_2-sort-trim_1.fq.gz,${INDIR}F1_PS_lib435254_7093_3/F1_PS_lib435254_7093_3-sort-trim_1.fq.gz \
--right ${INDIR}F1_EMIII_lib435236_7094_4/F1_EMIII_lib435236_7094_4-sort-trim_2.fq.gz,${INDIR}F1_EMII_lib435231_7093_2/F1_EMII_lib435231_7093_2-sort-trim_2.fq.gz,${INDIR}F1_EMI_lib435225_7093_2/F1_EMI_lib435225_7093_2-sort-trim_2.fq.gz,${INDIR}F1_ESII_lib435246_7093_3/F1_ESII_lib435246_7093_3-sort-trim_2.fq.gz,${INDIR}F1_ESI_lib435241_7093_2/F1_ESI_lib435241_7093_2-sort-trim_2.fq.gz,${INDIR}F1_PS_lib435254_7093_3/F1_PS_lib435254_7093_3-sort-trim_2.fq.gz \
--CPU 65 --max_memory 600G --output /home/gala0002/proj/proj_dir/2.0_trinity_out_F1/


##### ALL

nice -5 ${Trinity} --seqType fq --left ${INDIR}F1_EMIII_lib435236_7094_4-sort-trim_1.fq.gz --right ${INDIR}F1_EMIII_lib435236_7094_4-sort-trim_2.fq.gz --CPU 70 --max_memory 600G --output /home/gala0002/proj/proj_dir/2.0_trinity_out/


COMM
