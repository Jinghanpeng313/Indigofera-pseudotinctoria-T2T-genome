# Indigofera pseudotinctoria Haplotype-Resolved Genome Analyses Pipeline
### Analysis Pipelines and Plotting Codes

This repository contains the analysis pipelines and visualization scripts used in the study:

**“A haplotype-resolved T2T genome assembly of *Indigofera pseudotinctoria* reveals the genetic basis of flavonoid biosynthesis in Chinese Indigo.”**

The workflows cover **genome assembly, genome annotation, haplotype-resolved comparison, allele-specific expression (ASE), phylogenomics and whole-genome duplication (WGD), transcriptomic–metabolomic integration**, and **figure plotting** for the manuscript.


---

## Repository Structure

The repository contains seven major modules corresponding to the main analyses performed in this study:

```text
Indigofera-pseudotinctoria-T2T-genome/
├── 1.Genome_assembly/                      # T2T genome assembly pipeline
├── 2.Genome_annotation/                    # Structural and functional annotation
├── 3.Haplotype-resolved_genome_comparison/ # Comparative genomic analyses
├── 4.ASE_analysis/                         # Allele-specific expression analyses
├── 5.Phylogenomics_and_WGD/                # Phylogenomics and whole-genome duplication analyses
├── 6.Transcriptomic_metabolomic_analysis/  # Integrated transcriptome–metabolome analyses
└── 7.Plot_code/                            # Scripts for figure generation
```


### Module Description

- **`1.Genome_assembly/`**  
  Genome assembly workflows, including long-read assembly, Hi-C scaffolding, and assembly quality evaluation.

- **`2.Genome_annotation/`**  
  Structural and functional annotation of the haplotype-resolved genome assemblies (e.g. Hap1 and Hap2).

- **`3.Haplotype-resolved_genome_comparison/`**  
  Comparative analyses between haplotypes, including structural variations (SVs), transposable elements (TEs), and genome-wide divergence.

- **`4.ASE_analysis/`**  
  Allele-specific expression (ASE) analyses based on transcriptomic data.

- **`5.Phylogenomics_and_WGD/`**  
  Phylogenomic reconstruction, gene family evolution, and whole-genome duplication (WGD) analyses.

- **`6.transcriptomic–metabolomic_analysis/`**  
  Integrated transcriptomic and metabolomic analyses, with a focus on flavonoid biosynthesis pathways.

- **`7.plot_code/`**  
  R and Python scripts used to generate figures and supplementary plots for the manuscript.

> Each directory may contain multiple scripts. Please read the comments within individual scripts for detailed usage instructions.

---

## Computational Environment

The analyses were conducted using a combination of **Shell, Python, and R** scripts.

### Recommended Environment
- Python ≥ 3.8
- R ≥ 4.1

---

