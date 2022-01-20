## Prepare genomic regions

### Get chromosome sizes

<pre>samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
cut -f1,2 Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai|sort -t $'\t' -k1,1 > human_chromosomes.txt</pre>

### Get gene regions

<pre>awk -F '\t' '($3=="gene") {printf("%s\t%d\t%s\n",$1,int($4)-1,$5);}' Homo_sapiens.GRCh38.101.gtf | sort -t $'\t' -k1,1 -k2,2n | bedtools merge > human_genes.bed</pre>

### Get intergenic regions

<pre>bedtools complement -i human_genes.bed -g human_chromosomes.txt > human_intergenic.bed</pre>

### Get exonic regions (merged)

<pre>awk -F '\t' '($3=="exon") {printf("%s\t%d\t%s\n",$1,int($4)-1,$5);}' Homo_sapiens.GRCh38.101.gtf | sort -t $'\t' -k1,1 -k2,2n | bedtools merge > human_exons.bed</pre>

### Get intronic regions

<pre>bedtools complement -i <(cat human_intergenic.bed human_exons.bed | sort -t $'\t' -k1,1 -k2,2n) -g human_chromosomes.txt > human_introns.bed</pre>
