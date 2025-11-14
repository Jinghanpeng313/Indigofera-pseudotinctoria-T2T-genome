#!/usr/bin/bash
########1.Identify synteny genes from two genomes ---- gene level identification
python -m jcvi.formats.gff bed --type=mRNA --key=ID hap1.gff -o hap1.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID hap1.gff -o hap2.bed
python -m jcvi.compara.catalog ortholog hap1 hap2 --no_strip_names --cscore=0.99 --full
########2.Idenity Syteny region via anchorwave. Running Anchorwave
gmap_build --dir=./ --genomedb=hap1 ../hap1.fa
gmap_build --dir=./ --genomedb=hap2 ../hap2.fa
gmap -t 40 -A -f samse -d hap1 -D ./ ../hap1.cds.fa > hap1.gmap.self.sam
gmap -t 40 -A -f samse -d hap2 -D ./ ../hap1.cds.fa > hap2.vs.hap1.gmap.sami

anchorwave proali -t 40 \
    -i ../hap1.gff \
    -as ../hap1.cds.fa \
    -r ../hap1.fa \
    -a hap2.vs.hap1.gmap.sam \
    -ar hap1.gmap.self.sam \
    -s ../hap2.fa \
    -R 2 \
    -Q 2 \
    -n hap1_vs_hap2.gmap.anchor \
    -o hap1_vs_hap2.gmap.maf \
    -f hap1_vs_hap2.gmap.f.maf \
    -m 0

####3.Run diamond for calculate sequence similarity
diamond makedb --in hap1.pep -d hap1
diamond blastp -d hap1 -q hap2.pep -o hap1_vs_hap2.ultrasens.blast -p 20 --evalue 1e-5

#######4.Transform Anchorwave output to used format
python anchor.anchor2coord.py hap1_vs_Hhap2.gmap.anchor hap1_vs_hap2.gmap.anchor.coords
sort -k10,10 -k11,11 -k1,1n -k2,2n -k3,3n -k4,4n hap1_vs_hap2.gmap.anchor.coords > hap1_vs_hap2.gmap.maf.coords.sorted

#######5.Transferring Anchorwave coords file to sqlite3 database and blast file
sqlite3 hap1_vs_hap2.maf.db
create table dat(qstart INTEGER, qend INTEGER, sstart INTEGER, send INTEGER, qlens INTEGER, slens INTEGER,identity INTEGER, qstand text , sstand text , qchr text, schr text);
.separator "\t"
.import hap1_vs_hap2.gmap.maf.coords.sorted
.quit

#######6.Transferring blast file to sqlite3 database
sqlite3 hap1_vs_hap2.ultrasens.blast.db
create table dat(g1 text, g2 text, simarity INTEGER, length INTEGER, mismatch INTEGER, gapopen INTEGER,qstart INTEGER, qend INTEGER , sstart INTEGER , send INTEGER, evalue INTEGER, bitscore INTEGER);
.separator "\t"
.import hap1_vs_hap2.ultrasens.blast dat
.quit

#######7. Identification potential alleles
python3 Identify_alleles.py hap1.hap2.1x1.lifted.anchors coords.maf.db hap1_hap2.bed hap1_vs_hap2.ultrasens.blast.db hap1-hap2.alleles.Anchorwave.txt
