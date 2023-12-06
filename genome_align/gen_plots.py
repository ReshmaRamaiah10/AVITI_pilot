import os
import pegasus as pg
import pegasusio as io
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from subprocess import check_call
from scipy.stats import pearsonr
from adjustText import adjust_text


out_folder = "figure"
os.makedirs(out_folder, exist_ok=True)

gex_dataset = {
    "Ch7": {
        "illumina": "data/Ch7/illumina/raw_feature_bc_matrix.h5",
        "aviti": "data/Ch7/aviti/raw_feature_bc_matrix.h5",
    },
    "Ch10": {
        "illumina": "data/Ch10/illumina/raw_feature_bc_matrix.h5",
        "aviti": "data/Ch10/aviti/raw_feature_bc_matrix.h5",
    },
}

def gen_gex_plots(sample, chemistry):
    print(f"=============\nConsider sample {sample}:")
    data_il = pg.read_input(gex_dataset[sample]['illumina'])
    data_av = pg.read_input(gex_dataset[sample]['aviti'])

    # QC
    percent_mito = 10 if chemistry == 'fiveprime' else 20
    print(f"Use percent_mito = {percent_mito} % threshold.")
    pg.qc_metrics(data_il, min_genes=500, mito_prefix='MT-', percent_mito=percent_mito)
    pg.filter_data(data_il)

    pg.qc_metrics(data_av, min_genes=500, mito_prefix='MT-', percent_mito=percent_mito)
    pg.filter_data(data_av)

    # Plot Fig 2b
    cells_common = np.intersect1d(data_il.obs_names, data_av.obs_names)
    n_cells_common = cells_common.size
    n_cells_il_only = data_il.shape[0] - n_cells_common
    n_cells_av_only = data_av.shape[0] - n_cells_common

    df_fig2b = pd.DataFrame({'set': ['Both', 'AVITI only', 'Illumina only'], 'size': [n_cells_common, n_cells_av_only, n_cells_il_only]})

    cell_bar = sns.barplot(data=df_fig2b, x='set', y='size')
    for i in cell_bar.containers:
        cell_bar.bar_label(i,)
    cell_bar.set_xlabel("")
    cell_bar.set_ylabel("Number of cells")
    plt.tight_layout()
    plt.savefig(f"{out_folder}/{sample}_fig2b.png", dpi=500)
    plt.close()
    print("Figure 2b is generated!")

    # Plot Fig 2c and 2d
    data_il.obs["n_genes"] = data_il.X.getnnz(axis=1)
    data_av.obs["n_genes"] = data_av.X.getnnz(axis=1)

    data_il.obs["n_counts"] = data_il.X.sum(axis=1).A1
    data_av.obs["n_counts"] = data_av.X.sum(axis=1).A1

    df_qc_list = [
        pd.DataFrame({'n_genes': data_il.obs["n_genes"].values, 'n_counts': data_il.obs["n_counts"].values, 'technology': ["Illumina"] * data_il.shape[0]}),
        pd.DataFrame({'n_genes': data_av.obs["n_genes"].values, 'n_counts': data_av.obs["n_counts"].values, 'technology': ["AVITI"] * data_av.shape[0]}),
    ]
    df_qc = pd.concat(df_qc_list)

    median_n_genes = [data_il.obs["n_genes"].median(), data_av.obs["n_genes"].median()]
    median_n_counts = [data_il.obs["n_counts"].median(), data_av.obs["n_counts"].median()]

    gene_violin = sns.violinplot(data=df_qc, x='technology', y='n_genes')
    gene_violin.text(0+0.05, 3000, median_n_genes[0])
    gene_violin.text(1+0.05, 3000, median_n_genes[1])
    gene_violin.set_ylabel("Number of genes")
    gene_violin.set_xlabel("")
    plt.savefig(f"{out_folder}/{sample}_fig2c.png", dpi=500)
    plt.close()
    print("Figure 2c is generated!")

    umi_violin = sns.violinplot(data=df_qc, x='technology', y='n_counts')
    umi_violin.text(0+0.05, 10000, median_n_counts[0])
    umi_violin.text(1+0.05, 10000, median_n_counts[1])
    umi_violin.set_ylabel("Number of UMIs")
    umi_violin.set_xlabel("")
    plt.savefig(f"{out_folder}/{sample}_fig2d.png", dpi=500)
    plt.close()
    print("Figure 2d is generated!")

    # Calculate Pearson's Correlation on genes
    raw_expr_il = data_il.X.sum(axis=0).A1
    raw_expr_av = data_av.X.sum(axis=0).A1
    corr = pearsonr(raw_expr_il, raw_expr_av)
    print(corr)

    # Generate Fig 2e
    df_fig2e = pd.DataFrame({'Illumina': raw_expr_il, 'AVITI': raw_expr_av, 'gene': data_il.var_names})
    df_fig2e['Illumina'] = 1e6 * df_fig2e['Illumina'] / df_fig2e['Illumina'].sum()
    df_fig2e['AVITI'] = 1e6 * df_fig2e['AVITI'] / df_fig2e['AVITI'].sum()
    df_fig2e["Rat"] = np.abs(np.log(df_fig2e["Illumina"] + 10) - np.log(df_fig2e["AVITI"] + 10))
    df_fig2e["Illumina"] = np.log(df_fig2e["Illumina"] + 1)
    df_fig2e["AVITI"] = np.log(df_fig2e["AVITI"] + 1)
    df_fig2e["Diff"] = df_fig2e["Rat"] > np.log(2)
    df_fig2e.to_csv(f"{out_folder}/{sample}_fig2e_table.csv", index=False)

    check_call(['Rscript', 'plot_fig2e.R', f"{out_folder}/{sample}_fig2e_table.csv", f"{out_folder}/{sample}_fig2e.png", str(round(corr[0], 2))])

    print("=============")

if __name__ == "__main__":
    gen_gex_plots('Ch7', 'threeprime')
    gen_gex_plots('Ch10', 'fiveprime')
