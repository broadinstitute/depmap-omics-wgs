#!/bin/zsh

wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz

bedtools intersect -v -a hg38.1000.bed -b hg38-blacklist.v2.bed.gz \
        > hg38.1000.bl_excluded.2.bed
bedtools sort -i hg38.1000.bl_excluded.2.bed -g hg38.chrom.sizes \
        > hg38.1000.sorted.bed

# GC wig
bedtools nuc -fi Homo_sapiens_assembly38.fasta -bed hg38.1000.sorted.bed \
        > hg38.1000.gc.txt

cut -f1,2,3,5 hg38.1000.gc.txt \
        | sed '1d' \
                > hg38.1000.gc.bedgraph

bedGraphToBigWig hg38.1000.gc.bedgraph hg38.chrom.sizes hg38.1000.gc.bw
bigWigToWig hg38.1000.gc.bw hg38.1000.gc.wig

sed 's/#bedGraph section /fixedStep chrom=/g' hg38.1000.gc.wig \
        | sed 's/:/ start=/g' \
        | sed 's/\-[0-9]*/ step=1000 span=1000/g' \
        | awk '{if(index($1,"chr") == 1) {print $4 } else print $0}' \
                > hg38.1000.gc.formatted.wig

# Mappability wig
bwtool/bwtool summary -fill=0 hg38.1000.sorted.bed k100.Umap.MultiTrackMappability.bw stdout \
        | awk '{ print $1"\t"$2"\t"$3"\t"$6 }' \
                > hg38.1000.with_mean_mapscore.bed

bedGraphToBigWig hg38.1000.with_mean_mapscore.bed hg38.chrom.sizes hg38.1000.mapscore.bw

bigWigToWig hg38.1000.mapscore.bw hg38.1000.mapscore.wig

sed 's/#bedGraph section /fixedStep chrom=/g' hg38.1000.mapscore.wig \
        | sed 's/:/ start=/g' \
        | sed 's/\-[0-9]*/ step=1000 span=1000/g' \
        | awk '{if(index($1,"chr") == 1) {print $4 } else print $0}' \
                > hg38.1000.mapscore.formatted.wig
