#!/bin/bash

#########################################
# curate the TSS SAF
#########################################
wget -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
gtftools.py -t gencode.v29.annotation.gtf.tss.bed -w 1000 gencode.v29.annotation.gtf
awk '{OFS="\t"} {print $6,$1,$2,$3,$4}' gencode.v29.annotation.gtf.tss.bed > gencode.v29.annotation.gtf.tss.bed.saf

#########################################
# download files of interest
#########################################
grep '1 hour' metadata.tsv | grep bam | grep '100 nM' | grep ChIP-seq | grep GRCh38 | cut -f19,43 > sample_list.tsv

grep  -v 'dex' metadata.tsv | grep bam | grep -v 'eth' | grep GRCh38 | grep ChIP-seq | cut -f19,43 > sample_list2.tsv

cut -f1 sample_list.tsv | grep -wFf - sample_list2.tsv > sample_list2.tsv.tmp

cat sample_list2.tsv.tmp >> sample_list.tsv && rm sample_list2.tsv.tmp sample_list2.tsv

wget -N -i $(cut -f2 sample_list.tsv)


#atac
echo "https://www.encodeproject.org/files/ENCFF597SLV/@@download/ENCFF597SLV.bam
https://www.encodeproject.org/files/ENCFF248HDS/@@download/ENCFF248HDS.bam
https://www.encodeproject.org/files/ENCFF978DQZ/@@download/ENCFF978DQZ.bam
https://www.encodeproject.org/files/ENCFF020COS/@@download/ENCFF020COS.bam
https://www.encodeproject.org/files/ENCFF758ORC/@@download/ENCFF758ORC.bam
https://www.encodeproject.org/files/ENCFF809EKV/@@download/ENCFF809EKV.bam" > atac_samples.tsv

wget -N -i atac_samples.tsv

#########################################
# extract TSS counts
#########################################
for BAM in *bam ; do
  CNT=$BAM.cnt
  if [ ! -r $CNT ] ; then
    ../app/subread-1.6.4-source/bin/featureCounts  -a gencode.v29.annotation.gtf.tss.bed.saf -F SAF -T $(nproc) -o $BAM.cnt $BAM
  fi
done

