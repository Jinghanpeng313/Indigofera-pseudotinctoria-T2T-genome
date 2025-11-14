#!/usr/bin/env bash

################################################
# Step 1: Run Infernal cmscan to find Rfam families
################################################

# Define the path to the Rfam database
rfamdb="${database}/Rfam/latest/"  # Set the Rfam database path

# Get the genome size (in megabases) from the contig file, using esl-seqstat
# Z represents the genome size in megabases, we multiply by 2 and divide by 10^6 for conversion
Z=$(esl-seqstat --dna contig.fa | grep "Total " | tr -s [:space:] | sed 's/: /\t/g' | cut -f 2)
Zvalue=$(awk -v i="${Z}" 'BEGIN {print i*2/1000000}')

# Run cmscan to search for Rfam families
cmscan -Z "${Zvalue}" --rfam --nohmmonly --cut_ga --cpu 10 \
  -o contig.Rfam.txt --verbose --fmt 2 \
  --tblout contig.Rfam.tbl \
  --clanin "${rfamdb}/Rfam.clanin" \
  "${rfamdb}/Rfam.cm" contig.fa

# Process the output to remove overlaps
grep -v " = " contig.Rfam.tbl > contig.Rfam.deoverlapped.tbl  # Remove overlap lines

# Convert the cmscan output to GFF format
perl convert_cmscan_to_gff3.pl --cmscan --fmt2 --all --source Rfam \
  contig.Rfam.deoverlapped.tbl > contig.Rfam.deoverlapped.gff


###################################
# Step 2: rRNA Prediction with RNAmmer
###################################

# Running RNAmmer for rRNA prediction
# tsu = 5S or 5.8S rRNA, ssu = 16S or 18S rRNA, lsu = 23S or 28S rRNA
rnammer -S euk -m lsu,ssu,tsu -xml hap1.rRNA.rnammer.xml \
  -gff hap1.rRNA.rnammer.tbl -h hap1.rRNA.rnammer.hmmreport < hap1.fa

# The -S option specifies the organism type: arc = archaea, bac = bacteria, euk = eukaryotes.
# The -xml option specifies XML format output; -gff specifies GFF format output; -h outputs an HTML report.

# Convert RNAmmer results to GFF format
perl Convert_RNAmmer_to_gff3.pl --input contig.rRNA.rnammer.tbl > contig.rRNA.rnammer.gff


###################################
# Step 3: tRNA Prediction with tRNAscan-SE
###################################

# Running tRNAscan-SE to predict tRNAs
tRNAscan-SE -E -G -o contig.tRNA.out -f contig.tRNA.ss -m contig.tRNA.stats contig.fa
# -E option for eukaryotes; -G for genomic sequence; -o for output file; -f for secondary structure file; -m for statistics file

# Convert tRNAscan-SE output to GFF format
perl convert_tRNAScanSE_to_gff3.pl --input contig.tRNA.out > contig.tRNA.out.gff


###################################
# Step 4: Merge and organize data
###################################

# Filter out rRNA and tRNA from the Rfam GFF file
awk '$3 !~ "rRNA" && $3 !~ "tRNA" {print}' contig.Rfam.deoverlapped.gff > contig.Rfam.NorRNAtRNA.gff

# Combine tRNA, rRNA, and other ncRNA annotations into one GFF file
cat contig.tRNA.out.gff contig.rRNA.rnammer.gff contig.Rfam.NorRNAtRNA.gff > contig.ncRNA.gff

# Convert the combined GFF to a standardized GFF format and sort it
agat_convert_sp_gxf2gxf.pl -g contig.ncRNA.gff -o contig.ncRNA.final.gff

# Generate rRNA regions for protein gene filtering
awk 'BEGIN{OFS="\t"} $3=="rRNA" {print $1, $4, $5}' contig.ncRNA.final.gff > contig.rRNA.bed


###################################
# Step 5: LiftOver contig GFF to chromosome-level GFF
###################################

# Use LiftOver to convert contig-level GFF to chromosome-level GFF, if genome chain is provided
liftOver -gff contig.ncRNA.final.gff "${datadir}/genome.chain" genome.ncRNA.final.gff unmapped

# Perform ncRNA statistics analysis
perl Rfam.stat.v2.0.pl -in genome.ncRNA.final.gff -o genome.ncRNA.stat.xls
