# Note: bustools (current branch) is fine to use; For kallisto, must install the devel branch
# Note: These FASTQ files are technically Smart-seq3 data but we'll also pretend that they're bulk data when we analyze bulk
# Note: In the demultiplexed examples, we pretend data.R*.fastq.gz and data2.R*.fastq.gz are two completely different samples (in the non-demuxed examples; all the files are from a single experiment)
# Note: For non-demuxed data, we could supply our own whitelist or (like in this notebook) use bustools whitelist; for demuxed data, we'll need to completely skip these steps
# Note: bustools count must use -m for TCC (transcript-level) counting and must use --cm when there are no UMIs (like bulk data)
# Note: To get transcript-level expression via quant-tcc, we need to either a) Supply -l and -s ourselves (mean and SD of fragment length distribution), which is mandatory if reads are single-end; b) Supply a flens.txt file via -f (flens.txt is automatically generated in the output directory by kallisto for paired-end reads), or c) Do neither (i.e. don't do -f and don't do -l or -s either; in which case, all transcripts are assumed to have the exact same length when quant-tcc is run)
# Note: Should make getting transcript-level expression an optional option (because the algorithm currently takes a very long time)
# Note: --fr-stranded, --rf-stranded, and --unstranded are options that can be supplied to the "kallisto bus" command to specifically specify strand specificity (I didn't include them in the code here but should be options that we add to "kb count" for when it calls "kallisto bus")

### Download human index ###

wget --continue ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget --continue ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz Homo_sapiens.GRCh38.101.gtf.gz

### Download and unzip data files

wget --continue https://github.com/Yenaled/uploads/raw/main/data.I1.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data.I2.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data.R1.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data.R2.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data2.I1.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data2.I2.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data2.R1.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/data2.R2.fastq.gz
wget --continue https://github.com/Yenaled/uploads/raw/main/batch.txt
wget --continue https://github.com/Yenaled/uploads/raw/main/batch_single.txt


### Bulk (non-demuxed; paired-end) ###

# Remove existing output file
rm -rf output/
# Pseudoalignment, sort, whitelist, and correct
kallisto bus -x bulk -o output/ -i index.idx --paired data.I1.fastq.gz data.I2.fastq.gz data.R1.fastq.gz data.R2.fastq.gz data2.I1.fastq.gz data2.I2.fastq.gz data2.R1.fastq.gz data2.R2.fastq.gz
bustools sort -o output/output.s.bus output/output.bus
bustools whitelist -o output/whitelist.txt output/output.s.bus
bustools correct -o output/output.s.c.bus -w output/whitelist.txt output/output.s.bus
bustools sort -o output/output.unfiltered.bus output/output.s.c.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory)
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -i output/index.saved -f output/flens.txt -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx
# It's nice to have the barcodes in the quant/ directory too
cp output/counts_tcc/output.barcodes.txt output/quant/

### Bulk (non-demuxed; single-end) ###
# Same as above except cut down fastq files and use -l 200 -s 20 instead of -f flens.txt

# Remove existing output file
rm -rf output/
# Pseudoalignment, sort, whitelist, and correct
kallisto bus -x bulk -o output/ -i index.idx data.I1.fastq.gz data.I2.fastq.gz data.R1.fastq.gz data2.I1.fastq.gz data2.I2.fastq.gz data2.R1.fastq.gz
bustools sort -o output/output.s.bus output/output.bus
bustools whitelist -o output/whitelist.txt output/output.s.bus
bustools correct -o output/output.s.c.bus -w output/whitelist.txt output/output.s.bus
bustools sort -o output/output.unfiltered.bus output/output.s.c.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory); Note: here, for single-end, quant-tcc uses -l and -s instead of flens.txt
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -l 200 -s 20 -i output/index.saved -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx


### Bulk (demultiplexed; paired-end; files supplied on command-line) ###

# Remove existing output file
rm -rf output/
# Pseudoalignment, sort (no whitelist/barcode stuff because already demuxed)
kallisto bus -o output/ -i index.idx --paired data.R1.fastq.gz data.R2.fastq.gz data2.R1.fastq.gz data2.R2.fastq.gz
bustools sort -o output/output.unfiltered.bus output/output.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory)
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -i output/index.saved -f output/flens.txt -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx
# It's nice to have the barcodes in the quant/ directory too
cp output/counts_tcc/output.barcodes.txt output/quant/

### Bulk (demultiplexed; single-end;  files supplied on command-line) ###

# Remove existing output file
rm -rf output/
# Pseudoalignment, sort (no whitelist/barcode stuff because already demuxed)
kallisto bus -o output/ -i index.idx data.R1.fastq.gz data2.R1.fastq.gz
bustools sort -o output/output.unfiltered.bus output/output.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory)
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -l 200 -s 20 -i output/index.saved -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx
# It's nice to have the barcodes in the quant/ directory too
cp output/counts_tcc/output.barcodes.txt output/quant/

###  Bulk (demultiplexed; paired-end; files supplied in batch.txt file) ###

# Remove existing output file
rm -rf output/
# Pseudoalignment, sort (no whitelist/barcode stuff because already demuxed); Note: specifying --paired to "kallisto bus" is not required (though you can do so if you like) since kallisto can already infer single-/paired-end from batch file
kallisto bus -o output/ -i index.idx --batch=batch.txt
bustools sort -o output/output.unfiltered.bus output/output.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory)
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -i output/index.saved -f output/flens.txt -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx
# It's nice to have the barcodes in the quant/ directory too
cp output/counts_tcc/output.barcodes.txt output/quant/

###   Bulk (demultiplexed; single-end; files supplied in batch_single.txt file) ###

rm -rf output/
# Pseudoalignment, sort (no whitelist/barcode stuff because already demuxed)
kallisto bus -o output/ -i index.idx --batch=batch_single.txt
bustools sort -o output/output.unfiltered.bus output/output.bus
# Count (gene-level; put into counts_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
# Count (transcript-level; put into quant/ directory)
bustools count --cm -m -o output/counts_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output.unfiltered.bus
kallisto quant-tcc -o output/quant/ -l 200 -s 20 -i output/index.saved -e output/counts_tcc/output.ec.txt -g t2g.txt output/counts_tcc/output.mtx
# It's nice to have the barcodes in the quant/ directory too
cp output/counts_tcc/output.barcodes.txt output/quant/


### Smart-seq3 ###

rm -rf output/
# Pseudoalignment, sort, whitelist, and correct
kallisto bus -x smartseq3 -o output/ -i index.idx --paired data.I1.fastq.gz data.I2.fastq.gz data.R1.fastq.gz data.R2.fastq.gz data2.I1.fastq.gz data2.I2.fastq.gz data2.R1.fastq.gz data2.R2.fastq.gz
bustools sort -o output/output.s.bus output/output.bus
bustools whitelist -o output/whitelist.txt output/output.s.bus
bustools correct -o output/output.s.c.bus -w output/whitelist.txt output/output.s.bus
bustools sort -o output/output.unfiltered.bus output/output.s.c.bus
# Split into UMI and BULK(BULK records are identified by 32 T's in the UMI field) 
echo "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" > output/capture_nonUMI.txt
bustools capture -o output/output_internal.unfiltered.bus -c output/capture_nonUMI.txt -u output/output.unfiltered.bus
bustools capture -o output/output_umi.unfiltered.bus -c output/capture_nonUMI.txt -u -x output/output.unfiltered.bus
# Count UMI (gene-level; put into counts_umi_unfiltered/ directory); Since this is UMI, no need for --cm (note: smart-seq3 UMIs are short so need --umi-gene to handle them properly)
bustools count --umi-gene --genecounts -o output/counts_umi_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output_umi.unfiltered.bus
# Count BULK (gene-level; put into counts_internal_unfiltered/ directory)
bustools count --cm --genecounts -o output/counts_internal_unfiltered/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output_internal.unfiltered.bus
# Count UMI (transcript-level; put into quant_umi/ directory); Since this is UMI, no need for --cm (still always need -m when counting  TCCs though!) and no need for -f/-l/-s for quant-tcc
bustools count --umi-gene -m -o output/counts_umi_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output_umi.unfiltered.bus
kallisto quant-tcc -o output/quant_umi/ -i output/index.saved -e output/counts_umi_tcc/output.ec.txt -g t2g.txt output/counts_umi_tcc/output.mtx
cp output/counts_umi_tcc/output.barcodes.txt output/quant_umi/
# Count BULK (transcript-level; put into quant_internal/ directory); need -f flens.txt here
bustools count --cm -m -o output/counts_internal_tcc/ -g t2g.txt -t output/transcripts.txt -e output/matrix.ec output/output_internal.unfiltered.bus
kallisto quant-tcc -o output/quant_internal/ -i output/index.saved -f output/flens.txt -e output/counts_internal_tcc/output.ec.txt -g t2g.txt output/counts_internal_tcc/output.mtx
cp output/counts_internal_tcc/output.barcodes.txt output/quant_internal/
