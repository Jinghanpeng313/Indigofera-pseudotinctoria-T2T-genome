#!/usr/bin/env bash

#######################################################################
# Link Genome Files and Prepare Directories
#######################################################################

# Split the Hap1 and Hap2 genomes for multi-threaded annotation
perl split_fa_by_length.pl Hap1.fa Hap1.split
perl split_fa_by_length.pl Hap2.fa Hap2.split

# Create directory for BRAKER and navigate to it
mkdir -p 03.Braker
cd 03.Braker

# Set the number of threads for parallel processing
threads=10

# Define paths for RNA-seq data
fq1=$datadir/rnaseq_1.fq.gz
fq2=$datadir/rnaseq_2.fq.gz
genome=../Hap1.fa  # Use Hap1 for now, but will repeat for Hap2 in later steps

#######################################################################
# Step 1: HISAT2 and StringTie for Transcript Assembly
#######################################################################

# Build HISAT2 index for the Hap1 genome
hisat2-build Hap1 $genome

# Align RNA-seq reads (paired-end) to the Hap1 genome using HISAT2, then sort with SAMtools
hisat2 --dta --new-summary -p $threads -x $genome -1 $fq1 -2 $fq2 2>hisat2.log | samtools sort -@ 10 > rnaseq.bam
samtools index rnaseq.bam  # Index the BAM file for later visualization

# Assemble transcripts with StringTie
stringtie rnaseq.bam -p $threads -o rnaseq.gtf

#######################################################################
# Step 2: TransDecoder for ORF Prediction
#######################################################################

# Convert StringTie GTF to a fasta file containing predicted coding sequences (CDS)
cufflinks_gtf_genome_to_cdna_fasta.pl rnaseq.gtf $genome > transcripts.fasta

# Convert the GTF file to GFF3 format
cufflinks_gtf_to_alignment_gff3.pl rnaseq.gtf > transcripts.gff3

# Identify and predict ORFs in the transcript file
TransDecoder.LongOrfs -t transcripts.fasta
TransDecoder.Predict -t transcripts.fasta

# Map ORFs to genome and generate GFF3 file for downstream tools
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
sed 's#transdecoder#StringTie#' transcripts.fasta.transdecoder.genome.gff3 > stringtie.evm.format.gff3

#######################################################################
# Step 3: BRAKER Gene Annotation with Augustus and GeneMark-ET
#######################################################################

# Copy the GeneMark key file for use with GeneMark
cp /share/work/biosoft/GeneMarkS/gm_key_64 ~/.gm_key

# Run BRAKER for gene annotation on Hap1 genome using Augustus and GeneMark with RNA-seq data
braker.pl --species=$species \
   --genome=$genome \
   --softmasking \
   --bam=rnaseq.bam \
   --cores=$threads \
   --useexisting \
   --gff3 --skip_fixing_broken_genes

#######################################################################
# Step 4: PASA Assembly for Transcript Updates
#######################################################################

mkdir -p 04.PASA_assembly
cd 04.PASA_assembly

# Define paths to RNA-seq data and EST files
fq1=$datadir/rnaseq_1.fq.gz
fq2=$datadir/rnaseq_2.fq.gz
est=$datadir/est.fa  # Complete EST sequence for this species

# Trinity assembly using RNA-seq BAM file (genome-guided)
Trinity --genome_guided_bam ../03.Braker/rnaseq.bam \
        --genome_guided_max_intron 100000 \
        --max_memory 300G --CPU $threads \
        --output trinity_gg_out --full_cleanup

# Trinity de novo assembly (for higher gene density species)
Trinity --seqType fq --left $fq1 --right $fq2 --output trinity_tdn_out \
        --CPU $threads --full_cleanup --max_memory 300G

# Extract accession IDs from de novo assembly
accession_extractor.pl < trinity_tdn_out.Trinity.fasta > tdn.accs

# Combine de novo and genome-guided assemblies along with ESTs into a single transcript file
cat ../03.Braker/transcripts.fasta.transdecoder.cds trinity_gg_out.Trinity-GG.fasta trinity_tdn_out.Trinity.fasta $est > all_transcripts.fasta

#######################################################################
# Step 5: PASA Alignment and Filtering
#######################################################################

# Define the MySQL database for PASA
DATABASE=${species}_pasa
echo "
DATABASE=$DATABASE
" > alignAssembly.config

# Index and clean the transcript data
seqclean all_transcripts.fasta -c 10 -n 10000 -v $database/UniVec/UniVec

# Extract full-length transcript IDs
accession_extractor.pl < $est > FL_accs.txt

# Run PASA pipeline for transcript assembly comparison and integration
Launch_PASA_pipeline.pl -f FL_accs.txt \
    -c alignAssembly.config -C -r -R -g ../Hap1.fa \
    -T -t all_transcripts.fasta.clean -u all_transcripts.fasta \
    --CPU $threads --ALIGNERS blat,minimap2 --TRANSDECODER \
    --TDN tdn.accs --MAX_INTRON_LENGTH 1000000

# Second round of PASA pipeline for further updates
recent_update_gff_file=$(ls -t *.gene_structures_post_PASA_updates.*.gff3 | head -1)
Launch_PASA_pipeline.pl -f FL_accs.txt \
    -c alignAssembly.config -C -r -R -g ../Hap1.fa \
    -T -t all_transcripts.fasta.clean -u all_transcripts.fasta \
    --CPU $threads --ALIGNERS blat,minimap2 --TRANSDECODER \
    --TDN tdn.accs --MAX_INTRON_LENGTH 1000000

# Clean up and output the final GFF3 annotation
recent_update_gff_file=$(ls -t *.gene_structures_post_PASA_updates.*.gff3 | head -1)
cat $recent_update_gff_file | grep -v "^#" | sed "/^$/d" > genome.all.gff3

# Filter out ORFs shorter than 50 amino acids
orflen=100
agat_sp_filter_by_ORF_size.pl --gff genome.all.gff3 --size $orflen --test ">" --output genome.all${orflen}.gff3

# Filter incomplete genes
agat_sp_filter_incomplete_gene_coding_models.pl --gff genome.all${orflen}3_sup${orflen}.gff --fasta ../genome.fa --output genome.filtered_partial_gene.gff3

# Filter out non-functional genes based on annotations
ln -s genome.filtered_partial_gene.gff3 genome.filtered_func.gff3

# Keep the longest transcript for each gene
agat_sp_keep_longest_isoform.pl --gff genome.filtered_func.gff3 -o genome.final.longest_isoform.gff3

# Filter genes overlapping with rRNA regions
agat_sp_filter_record_by_coordinates.pl --gff genome.final.longest_isoform.gff3 -c ../02.ncRNA/genome.rRNA.bed -o genome.rRNA-filter

# Convert genome-level GFF to genome-level GFF using liftOver
liftOver -gff genome.rRNA-filter/remaining.gff3 $datadir/genome.chain genome.liftOver.gff unmapped

# Rename gene IDs for clarity
agat_sp_manage_IDs.pl --gff genome.liftOver.gff --prefix ${species} --ensembl --collective --type_dependent -o genome.final.gff3

# Generate protein sequences from the final GFF annotation
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl --gff3 genome.final.gff3 --fasta ../genome.fa --seqType cDNA > cdna.fa
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl --gff3 genome.final.gff3 --fasta ../genome.fa --seqType CDS > cds.fa
/share/work/biosoft/TransDecoder/latest/util/gff3_file_to_proteins.pl --gff3 genome.final.gff3 --fasta ../genome.fa --seqType prot > pep.fa

#######################################################################
# Repeat the Above Steps for Hap2 (You can copy and modify the steps)
#######################################################################

# Use Hap2.fa as the genome input for the same set of steps for Hap2
genome=../Hap2.fa

# The rest of the steps (HISAT2, StringTie, TransDecoder, BRAKER, PASA, etc.) can be repeated similarly for Hap2.fa

#######################################################################
# End of Script
#######################################################################
