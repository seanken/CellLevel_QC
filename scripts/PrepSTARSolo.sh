bam=$1 ##The input bam from STARSolo
gtf=$2 ##The GTF to use
outprefix=$3 ##The name of the output prefix to use

echo Make gtf for labelling

genes=${outprefix}.genes.gtf ##GTF with gene info
exons=${outprefix}.exons.gtf ##GTF with exon info
anti=${outprefix}.anti.gtf ##GTF with genes info but opposite strand


awk '{if($3=="gene"){print $0}}' $gtf > $genes
awk '{if($3=="exon"){print $0}}' $gtf > $exons
awk -F "\t" -v OFS="\t" '{if($7=="+"){$7="-"}else{$7="+"} print $0}' $genes > $anti

echo Label bam
bedtools tag -s -tag RE -i $bam -files $genes $exons $anti -labels N E A > ${outprefix}.label.bam
