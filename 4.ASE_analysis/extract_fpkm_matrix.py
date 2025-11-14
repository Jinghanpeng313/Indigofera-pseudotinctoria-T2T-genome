import os
import pandas as pd
import re
from glob import glob
from tqdm import tqdm

gtf_dir = "/home/jinghanpeng/data/MJ_analysis/06ASE/hap2/gtf_files/"
gtf_files = sorted(glob(os.path.join(gtf_dir, "*.transcripts.gtf")))

fpkm_dict = {}

gene_fpkm_pattern = re.compile(r'gene_id "([^"]+)";.*?FPKM "([^"]+)"')

for gtf in tqdm(gtf_files, desc="处理GTF文件"):
    sample = os.path.basename(gtf).split(".")[0]
    gene_fpkm = {}

    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" not in line:
                continue
            match = gene_fpkm_pattern.search(line)
            if match:
                gene_id, fpkm = match.groups()
                fpkm = float(fpkm)
                gene_fpkm[gene_id] = fpkm

    fpkm_dict[sample] = gene_fpkm


df = pd.DataFrame.from_dict(fpkm_dict, orient='index').fillna(0).T
df.index.name = "GeneID"


df.to_csv("all_treat_expre_filter_fpkm_for_pca.txt", sep="\t")
print("FPKM martix all_treat_expre_filter_fpkm_for_pca.txt")
