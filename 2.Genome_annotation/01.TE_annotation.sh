#!/usr/bin/env bash
# ---------- Load environment variables or read from ENV ----------
if [[ -f "./env.sh" ]]; then
  # env.sh：workdir, datadir, database, scriptsdir, threads, TpasesPROT, SPROT
  # export workdir="/path/to/work"
  # export datadir="/path/to/data"
  # export TpasesPROT="/db/Tpases020812"
  # export SPROT="/db/uniprot_sprot_clean.fasta"
  # export threads=16
  source ./env.sh
else
  log "env.sh not found. Falling back to environment variables."
fi

: "${workdir:?need workdir}"
: "${datadir:?need datadir}"
: "${database:?need database}"
: "${scriptsdir:?need scriptsdir}"
: "${threads:=10}"
: "${TpasesPROT:?need TpasesPROT transposase DB}"
: "${SPROT:?need SPROT (UniProt-SwissProt without transposases)}"


# ---------- Prepare directories and files ----------

REP_DIR="${workdir}/01.Repeat"
mkcd() {
  mkdir -p "$1" && cd "$1"
}

mkcd "${REP_DIR}"

ln -sf "${datadir}/Hap1.fasta" Hap1.fa
ln -sf "${datadir}/genome.fasta" genome.fa
# seqkit replace -p '(\S+)' -r 'ctg_{nr}' Hap1.fa > Hap1.simplified.fa

# ---------- Step 1: LTR prediction using ltrharvest, LTR_FINDER_parallel, and LTR_retriever ----------
log "Step 1: LTR prediction (ltrharvest, LTR_FINDER_parallel) and LTR_retriever"

mkcd "${REP_DIR}/01.ltrharvest"
ln -sf ../Hap1.fa .
gt suffixerator -db Hap1.fa -indexname Hap1.fa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index Hap1.fa -out Hap1.fa.out85 -outinner Hap1.fa.outinner85 -gff3 Hap1.fa.gff85 \
  -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 \
  -mintsd 5 -maxtsd 5 -vic 10 -similar 85 > Hap1.fa.harvest.scn

mkcd "${REP_DIR}/02.ltrfinder"
ln -sf ../Hap1.fa .
LTR_FINDER_parallel -seq Hap1.fa -threads "${threads}" -harvest_out -size 1000000 -time 300

mkcd "${REP_DIR}/03.LTR_retriever"
ln -sf ../Hap1.fa .
cat ../01.ltrharvest/Hap1.fa.harvest.scn ../02.ltrfinder/Hap1.fa.finder.combine.scn > Hap1.fa.rawLTR.scn
LTR_retriever -genome Hap1.fa -inharvest Hap1.fa.rawLTR.scn -threads "${threads}" -u 1.3e-8

# LTR insertion time plot
Rscript ltr_divergence_time_plot.r -i Hap1.fa.pass.list -p ltr_insert_time

# ---------- Step 2: Run RepeatModeler on masked Hap1s ----------
log "Step 2: RepeatModeler on masked Hap1s"

mkcd "${REP_DIR}/05.RepeatModeler"
ln -sf ../Hap1.fa 
cat  ../03.LTR_retriever/Hap1.fa.LTRlib.fa > MITE_LTR.lib

RepeatMasker -pa "${threads}" -no_is -lib MITE_LTR.lib -dir . Hap1.fa
BuildDatabase -name um_Hap1.db -engine rmblast Hap1.fa.masked
RepeatModeler -pa "${threads}" -database um_Hap1.db >& RepeatModeler.log

perl "CRL_Scripts1.0/repeatmodeler_parse.pl" --fastafile RM*/consensi.fa.classified \
  --unknowns repeatmodeler_unknowns.fasta --identities repeatmodeler_identities.fasta

# BLAST unknown sequences against TpasesPROT database
blastx -query repeatmodeler_unknowns.fasta -db "${TpasesPROT}" -evalue 1e-10 -num_descriptions 10 \
  -out modelerunknown_blast_results.txt -num_threads "${threads}"

perl "CRL_Scripts1.0/transposon_blast_parse.pl" \
  --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknowns.fasta

ln -sf unknown_elements.txt ModelerUnknown.lib
cat identified_elements.txt repeatmodeler_identities.fasta > ModelerKnown.lib
cat ModelerUnknown.lib ModelerKnown.lib > ModelerAll.lib

# ---------- Step 3: Homology-based masking ----------
log "Step 3: Homology-based masking (RepeatMasker lineage)"

mkcd "${REP_DIR}/06.RepeatMasker-homo"
ln -sf ../Hap1.fa .
RepeatMaskerLib="$(python - <<'PY'
import RepeatMasker as rm, os
print(getattr(rm, "LIB_DIR", "/usr/share/RepeatMasker/Libraries") + "/RepeatMaskerLib.h5")
PY
)"
famdb.py -i "${RepeatMaskerLib}" lineage --ancestors --descendants Sapindales || true
famdb.py -i "${RepeatMaskerLib}" families --format fasta_name --ancestors --curated --descendants \
  --include-class-in-name "Sapindales" > Homology.db

RepeatMasker -e rmblast -gff -xsmall -no_is -nolow -pa "${threads}" -norna -lib Homology.db Hap1.fa

# ---------- Step 4: Merge libraries and deduplicate ----------
log "Step 4: Merge libraries, ProtExcluder, deduplicate, DeepTE classify"

mkcd "${REP_DIR}/07.RepeatMasker"
ln -sf ../Hap1.fa .
ln -sf ../05.RepeatModeler/ModelerAll.lib .
ln -sf ../06.RepeatMasker-homo/Homology.db .

for lib in ModelerAll.lib MITE_LTR.lib Homology.db; do
  blastx -query "${lib}" -db "${SPROT}" -evalue 1e-10 -num_descriptions 10 -num_threads "${threads}" \
         -out "${lib}_blast_results.txt"
  perl "ProtExcluder1.2/ProtExcluder.pl" "${lib}_blast_results.txt" "${lib}"
  log "${lib} before=$(grep -c '>' ${lib} || echo 0) after=$(grep -c '>' ${lib}noProtFinal || echo 0)"
done

cat ModelerAll.libnoProtFinal Homology.dbnoProtFinal > allRepeats.lib
vsearch --cluster_fast allRepeats.lib --id 0.97 --centroids allRepeats.non-redundant.lib

# ---------- Step 5: Mask genome and generate RepeatLandscape ----------
log "Step 5: Mask genome with final library & landscape plot"

RepeatMasker -e rmblast -gff -xsmall -html -norna -source -no_is -pa "${threads}" -s \
             -lib allRepeats.final.lib Hap1.fa

Hap1_size=$(esl-seqstat --dna Hap1.fa | awk -F':' '/^Total /{gsub(/ /,"",$2);print $2}')
calcDivergenceFromAlign.pl -s Hap1.fa.divsum Hap1.fa.cat.gz
createRepeatLandscape.pl -div Hap1.fa.divsum -g "${Hap1_size}" > Hap1.fa.html

sed -n '1,/^$/p' Hap1.fa.divsum > Hap1.plot.identity
Rscript "kimura_lineplot.r" -i Hap1.plot.identity -p repeat_class_identity

sed -n '/Div/,$s/ $//p' Hap1.fa.divsum > Hap1.plot.divsum
Rscript "kimura_divergence_barplot.r" -i Hap1.plot.divsum -p repeat_class_divergence -g "${Hap1_size}"

perl "repeat_to_gff.pl" -in Hap1.fa.out -out Hap1.TE_repeat.gff
bedtools maskfasta -fi ../Hap1.fa -bed Hap1.TE_repeat.gff -fo Hap1.hardmasked.fa
bedtools maskfasta -soft -fi ../Hap1.fa -bed Hap1.TE_repeat.gff -fo Hap1.softmasked.fa

# ---------- Step 6: RepeatProteinMask ----------
log "Step 6: RepeatProteinMask on hardmasked Hap1s"

mkcd "${REP_DIR}/08.RepeatProteinMask"
ln -sf ../07.RepeatMasker/Hap1.hardmasked.fa Hap1.fa
rm -f cmds.sh
perl "split_fa_by_length.pl" Hap1.fa Hap1.split
num=$(ls Hap1.split/seq*fa | wc -l)
for i in $(seq 1 "${num}"); do
  echo "RepeatProteinMask -e rmblast -noLowSimple -pvalue 0.0001 Hap1.split/seq${i}.fa" >> cmds.sh
done

ParaFly -c cmds.sh -CPU "${threads}"
cat Hap1.split/*annot | sed '/^pValue.*/d' > Hap1.RepeatProteinMask.annot
perl "repeat_to_gff.pl" -in Hap1.RepeatProteinMask.annot -out Hap1.RepeatProteinMask.annot.gff
awk 'BEGIN{sum=0}{sum+=($5-$4+1)}END{print "RepeatProteinMask_length:",sum}' Hap1.RepeatProteinMask.annot.gff

# ---------- Step 7: TRF simple repeats and final repeat masking ----------
log "Step 7: TRF simple repeats & merge"

mkcd "${REP_DIR}/09.RepeatFinal"
cat ../08.RepeatProteinMask/Hap1.RepeatProteinMask.annot.gff \
    ../07.RepeatMasker/Hap1.TE_repeat.gff > Hap1.all.TE_repeat.gff

bedtools maskfasta -fi ../Hap1.fa -bed Hap1.all.TE_repeat.gff -fo Hap1.hardmasked.fa
bedtools maskfasta -soft -fi ../Hap1.fa -bed Hap1.all.TE_repeat.gff -fo Hap1.softmasked.fa

trf Hap1.hardmasked.fa 2 7 7 80 10 50 2000 -h -d
perl "repeat_to_gff.pl" \
  -in Hap1.hardmasked.fa.2.7.7.80.10.50.2000.dat -out Hap1.trf.gff
awk 'BEGIN{sum=0}{sum+=($5-$4+1)}END{print "TRF_length:",sum}' Hap1.trf.gff

cat Hap1.all.TE_repeat.gff Hap1.trf.gff > Hap1.all.repeat.gff

# ---------- Step 8: LiftOver (Optional) ----------
log "Step 8: LiftOver to chromosome-level (if chain file is provided)"

liftOver -gff Hap1.all.repeat.gff "${datadir}/genome.chain" genome.all.repeat.gff unmapped || true
liftOver -gff Hap1.all.TE_repeat.gff "${datadir}/genome.chain" genome.TE_repeat.gff unmapped || true

bedtools maskfasta -fi ../genome.fa -bed genome.TE_repeat.gff -fo genome.hardmasked.fa
bedtools maskfasta -soft -fi ../genome.fa -bed genome.TE_repeat.gff -fo genome.softmasked.fa

samtools faidx ../genome.fa
perl "buildRepeatSummary.pl" -genome ../genome.fa.fai -useAbsoluteGenomeSize genome.all.repeat.gff \
  > genome.repeat.summary.txt

log "All steps finished successfully."
