#!/usr/bin/env bash
# ----------------------------
# Configurable parameters
# ----------------------------
THREADS=48
HIFIASM_THREADS=64
WORKDIR=$(pwd)
PREFIX="sample"
HIFI_READS="hifi_reads.fq.gz"
ONT_ULTRA="ont_ultra_long.fq.gz"
ONT_RAW="ont_raw_reads.fastq.gz"
UL_FILE="ul.fq.gz"
HIC_R1="hic_R1.fq.gz"
HIC_R2="hic_R2.fq.gz"
ASSEMBLY_REF="${PREFIX}.asm.fa"   # placeholder for downstream steps

log() { echo -e "[`date '+%F %T'`] $*"; }

# ----------------------------
# Stage 1: Primary haplotype assembly with hifiasm
# ----------------------------
log "Stage 1: Running hifiasm for haplotype-resolved assembly (HiFi + optional ultra-long)"
# Note: adjust hifiasm flags to match your hifiasm version. --ul accepts a FASTA/FASTQ of ultra-long reads.
# If you have parental short-reads for trio-binning use --h1/--h2; for Hi-C based phasing newer hifiasm versions accept --hic
hifiasm -o ${PREFIX}_hifiasm -t ${HIFIASM_THREADS} \
  --ul ${UL_FILE} \
  ${HIFI_READS}

# hifiasm usually produces several outputs; find haplotig GFA/FASTA files
# Common names: <prefix>.bp.p_ctg.gfa ; <prefix>.bp.hap1.p_ctg.gfa ; <prefix>.bp.hap2.p_ctg.gfa
log "Stage 1: Extract haplotype contigs from GFA (if present)"
# Prefer gfatools if available; else fallback to awk extraction of S-lines (3rd column is sequence in many GFA versions)
if command -v gfatools >/dev/null 2>&1; then
  if [[ -f "${PREFIX}_hifiasm.bp.hap1.p_ctg.gfa" ]]; then
    gfatools gfa2fa ${PREFIX}_hifiasm.bp.hap1.p_ctg.gfa > haplotype1.raw.fa
  fi
  if [[ -f "${PREFIX}_hifiasm.bp.hap2.p_ctg.gfa" ]]; then
    gfatools gfa2fa ${PREFIX}_hifiasm.bp.hap2.p_ctg.gfa > haplotype2.raw.fa
  fi
else
  # awk fallback: works if S lines have sequence in 3rd column
  for f in ${PREFIX}_hifiasm.bp.hap1.p_ctg.gfa ${PREFIX}_hifiasm.bp.hap2.p_ctg.gfa; do
    if [[ -f "$f" ]]; then
      out="$(basename $f .gfa).fa"
      awk -F'\t' '/^S/{print ">"$2"\n"$3}' "$f" > "${out}"
    fi
  done
fi

# If hifiasm produced .p_ctg.fa directly, copy/rename
for x in ${PREFIX}_hifiasm*.p_ctg.fa ${PREFIX}_hifiasm*.hap*.fa 2>/dev/null; do
  [[ -f "$x" ]] || continue
  case "$x" in
    *hap1* ) cp -n "$x" haplotype1.raw.fa ;;
    *hap2* ) cp -n "$x" haplotype2.raw.fa ;;
    *p_ctg.fa) cp -n "$x" "${PREFIX}.primary_contigs.fa" ;;
  esac
done

# Basic check
for H in haplotype1.raw.fa haplotype2.raw.fa; do
  if [[ -f "$H" ]]; then
    log "Found: $H (size: $(du -h $H | cut -f1))"
  else
    log "Warning: $H not found — check hifiasm outputs."
  fi
done

# ----------------------------
# Stage 2: ONT ultra-long preliminary assembler (NextDenovo)
# ----------------------------
log "Stage 2: ONT ultra-long preliminary assembly with NextDenovo (as optional long-read draft)"
# You need a nextdenovo config file. We write a minimal template if not present.
if [[ ! -f nextdenovo.cfg ]]; then
  cat > nextdenovo.cfg <<'NDCFG'
[General]
job_type = local
job_prefix = nextDenovo
task = all
input_type = raw
read_type = ont
input_fofn = ./ont.fofn
workdir = nextdenovo_out
NDCFG
  # create fofn
  echo "${ONT_ULTRA}" > ont.fofn
fi

# run NextDenovo (change command if your nextDenovo binary name differs)
if command -v nextDenovo >/dev/null 2>&1; then
  nextDenovo ./nextdenovo.cfg || log "nextDenovo failed - check logs"
else
  log "nextDenovo not found in PATH; skipping ONT assembly step."
fi

# ----------------------------
# (Optional) Stage 3: Purge haplotigs / dedup (purge_dups) and scaffolding prep
# ----------------------------
log "Stage 3: (recommended) Purging redundants and preparing for scaffolding"
# Example: use purge_dups to remove allelic contigs before scaffolding (optional)
# purge_dups steps are dataset-specific and require read coverage histogram -> user should run separately as needed.

# ----------------------------
# Stage 4: Hi-C mapping & haplotype-aware validation/scaffolding
# ----------------------------
log "Stage 4: Map Hi-C reads to each haplotype and prepare filtered BAMs for HapHiC/3D-DNA/juicer"
# Build index for each haplotype assembly and map
for HAP in haplotype1.raw.fa haplotype2.raw.fa; do
  if [[ ! -f "$HAP" ]]; then
    log "Skipping Hi-C mapping for $HAP (file not present)"
    continue
  fi
  base=$(basename "$HAP" .raw.fa)
  # bwa index if not present
  if [[ ! -f "${HAP}.bwt" ]]; then
    bwa index "$HAP"
  fi

  # Map paired Hi-C reads (use bwa mem -5SP for Hi-C; check bwa version)
  bwa mem -5SP -t ${THREADS} "$HAP" "${HIC_R1}" "${HIC_R2}" \
    | samblaster -e -i --addMateTags \
    | samtools view -b -@ ${THREADS} -F 3340 - \
    | samtools sort -@ ${THREADS} -o ${base}.hic.sorted.bam -

  samtools index ${base}.hic.sorted.bam
  log "Generated ${base}.hic.sorted.bam"
done

# If you have HapHiC utilities, run haplotype-specific pipeline (paths may vary)
if [[ -x "HapHiC/haphic" ]]; then
  log "Running HapHiC pipeline (if available)"
  [[ -f haplotype1.raw.fa ]] && HapHiC/haphic pipeline haplotype1.raw.fa haplotype1.raw.fa.hic.sorted.bam 8 --threads ${THREADS} --processes ${THREADS} || true
  [[ -f haplotype2.raw.fa ]] && HapHiC/haphic pipeline haplotype2.raw.fa haplotype2.raw.fa.hic.sorted.bam 8 --threads ${THREADS} --processes ${THREADS} || true
else
  log "HapHiC not found or executable missing; skip HapHiC step. Consider using juicer/3D-DNA or SALSA for scaffolding."
fi

# ----------------------------
# Stage 5: ONT-based gap closure (TGS-GapCloser via Docker + alternative path)
# ----------------------------
log "Stage 5: Gap closing with TGS-GapCloser (Docker) and other ONT-based fillers"
for H in haplotype1.raw.fa haplotype2.raw.fa; do
  if [[ ! -f "$H" ]]; then continue; fi
  base=$(basename "$H" .raw.fa)
  out="${base}_gapfilled.fasta"

  # Using Docker image (example); adapt mount paths to your environment.
  if docker image inspect registry.cn-hangzhou.aliyuncs.com/kilianhett/tgsgapcloser:v3 >/dev/null 2>&1; then
    docker run --rm -v "${WORKDIR}":/data registry.cn-hangzhou.aliyuncs.com/kilianhett/tgsgapcloser:v3 \
      tgsgapcloser -l /data/${H} -s /data/${ONT_ULTRA} -o /data/${out} --minmap_arg "-x asm20" -t 32 -b 5000 || log "tgsgapcloser failed for ${H}"
  else
    log "Docker image for TGS-GapCloser not available locally; skipping docker step for ${H}"
  fi
done

# Quartet/NextDenovo/Hicanu combined gap filling (user provided scripts)
log "Stage 5b: Run Hicanu + Quartet GapFiller if available"
if command -v hicanu >/dev/null 2>&1; then
  mkdir -p hicanu_output
  # example hicanu invocation - adjust genomeSize and options
  hicanu -p ont_assembly -d ./hicanu_output genomeSize=0.66g -nanopore-raw ${ONT_RAW} \
    correctedErrorRate=0.5 batOptions="-trim-reads overlapEndTrimLimit=1000" || log "hicanu failed"
  cp ./hicanu_output/ont_assembly.contigs.fasta ./final_contigs.fasta || true
else
  log "hicanu not found; skip Hicanu step"
fi

# Example quartet.py GapFiller call (requires quartet.py in PATH)
if [[ -f "quartet.py" ]] && [[ -f "./hicanu_output/ont_assembly.contigs.fasta" ]]; then
  DRAFT_GENOME="./hicanu_output/ont_assembly.contigs.fasta"
  HICANU_CONTIG="./hicanu_output/ont_assembly.purged.fasta"
  NEXTDENOVO_CONTIG="./nextdenovo_out/assembly.fasta"
  python3 quartet.py GapFiller -d ${DRAFT_GENOME} -g ${HICANU_CONTIG} ${NEXTDENOVO_CONTIG} -f 10000 -l 2000 -i 60 -t 16 --minimapoption "-x asm20" --overwrite -p output_gapfilled || true
fi

# ----------------------------
# Stage 6: HiFi-based polishing (NextPolish or other polishers)
# ----------------------------
log "Stage 6: Polishing with NextPolish (HiFi reads) — ensure nextpolish.cfg exists and points to reads/indexes"
if [[ ! -f nextpolish.cfg ]]; then
  cat > nextpolish.cfg <<'NPC'
[General]
job_type = local
job_prefix = nextpolish
task = all
rewrite = yes

[sgs_options]
sgs_1 = hifi_reads.fq.gz
NPC
  log "Wrote template nextpolish.cfg — edit it to match your HiFi read files"
fi

for H in haplotype1.raw.fa haplotype2.raw.fa; do
  base=$(basename "$H" .raw.fa)
  gapfilled="${base}_gapfilled.fasta"
  polish_out="${base}.polished.fasta"

  # prefer gapfilled output if present, else use raw
  input_fa="$H"
  [[ -f "${gapfilled}" ]] && input_fa="${gapfilled}"

  if command -v next_polish >/dev/null 2>&1 || command -v nextpolish >/dev/null 2>&1; then
    # nextPolish sometimes called nextpolish or next_polish depending on install
    NEXTPOLISH_CMD="$(command -v nextpolish || command -v next_polish)"
    ${NEXTPOLISH_CMD} -c nextpolish.cfg -p ${THREADS} ${input_fa} -o ${polish_out} || log "nextpolish failed for ${input_fa}"
  else
    log "nextpolish not available; suggest using Racon/Medaka (for ONT) and Arrow/DeepVariant (for HiFi) as alternatives."
  fi
done

# ----------------------------
# Stage 7: Haplotype-specific QC: BUSCO, LTR analysis, basic stats
# ----------------------------
log "Stage 7: QC - BUSCO, LTR & assembly stats"
for H in haplotype1.polished.fasta haplotype2.polished.fasta; do
  if [[ ! -f "$H" ]]; then
    # try fallback names
    H_alt="${H/.polished.fasta/.polished.fasta}"
    [[ -f "$H_alt" ]] && H="$H_alt"
  fi
  if [[ -f "$H" ]]; then
    base=$(basename "$H" .polished.fasta)
    # BUSCO
    if command -v busco >/dev/null 2>&1; then
      busco -i ${H} -l embryophyta_odb10 -o ${base}_busco -m genome --cpu ${THREADS} || log "BUSCO failed on ${H}"
    else
      log "BUSCO not found; skipping"
    fi

    # LTR detection (example commands; adjust as needed)
    if command -v ltr_finder >/dev/null 2>&1 && command -v ltr_retriever >/dev/null 2>&1; then
      ltr_finder ${H} > ${base}.ltr_finder.out || log "ltr_finder failed on ${H}"
      ltr_retriever -genome ${H} -inharvest ${base}.ltr_finder.out -threads ${THREADS} -seq_len 10000000 || log "ltr_retriever failed on ${H}"
    else
      log "ltr_finder/ltr_retriever not available; skipping LTR analysis"
    fi

    # assembly stats
    if command -v seqkit >/dev/null 2>&1; then
      seqkit stats ${H} > ${base}.seqkit.stats.txt
    else
      awk '/^>/{if(N){print N; N=0}}; /^[^>]/ {N+=length($0)} END{if(N)print N}' ${H} > ${base}.lengths.txt || true
    fi
  else
    log "Polished file ${H} not found; skip QC"
  fi
done

log "Pipeline finished (or reached end of provided steps). Review logs and tool outputs for errors."
