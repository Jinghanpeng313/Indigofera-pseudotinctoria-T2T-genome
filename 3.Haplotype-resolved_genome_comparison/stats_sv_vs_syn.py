import pandas as pd
from scipy.stats import mannwhitneyu

# TPM matrix: rows = gene_id, columns = samples
tpm = pd.read_csv("tpm_matrix.tsv", sep="\t", index_col=0)

# Group gene lists
sv_ids  = pd.read_csv("hap1_SV.ids",  header=None)[0].tolist()
syn_ids = pd.read_csv("hap1_SYN_only.ids", header=None)[0].tolist()

# Keep only genes that exist in expression matrix
sv_ids  = [g for g in sv_ids  if g in tpm.index]
syn_ids = [g for g in syn_ids if g in tpm.index]

# Mean TPM per gene across samples
gene_mean = tpm.mean(axis=1)

sv_expr  = gene_mean.loc[sv_ids]
syn_expr = gene_mean.loc[syn_ids]

# Two-sided Mann–Whitney U test
u_stat, p_val = mannwhitneyu(sv_expr.values, syn_expr.values, alternative="two-sided")

print(f"SV genes (n={len(sv_expr)}), Syntenic genes (n={len(syn_expr)})")
print(f"Mann–Whitney U two-sided p-value: {p_val:.3e}")

# Optional: export data for plotting
out = pd.DataFrame({
    "gene_id": sv_expr.index.tolist() + syn_expr.index.tolist(),
    "group":   ["SV"]*len(sv_expr) + ["Syntenic"]*len(syn_expr),
    "mean_TPM": pd.concat([sv_expr, syn_expr]).values
})
out.to_csv("sv_vs_syn_meanTPM.tsv", sep="\t", index=False)
