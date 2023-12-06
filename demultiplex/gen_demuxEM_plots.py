import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pegasus as pg

from pegasusio import MultimodalData


pegasus_20 = [
    "#c5b0d5",
    "#ff7f0e",
    "#8c564b",
    "#ff9896",
    "#1f77b4",
    "#dbdb8d",
    "#e377c2",
    "#2ca02c",
    "#aec7e8",
    "#ffbb78",
    "#9edae5",
    "#98df8a",
    "#d62728",
    "#9467bd",
    "#c49c94",
    "#f7b6d2",
    "#bcbd22",
    "#17becf",
    "#ad494a",
    "#8c6d31",
]
# sample to donor data
assignment2donor = {
    "CMO301": "Donor 1",
    "CMO302": "Donor 1",
    "CMO303": "Donor 2",
    "CMO304": "Donor 2",
    "CMO305": "Donor 3",
    "CMO306": "Donor 3",
    "CMO307": "Donor 4",
    "CMO308": "Donor 4",
    "CMO309": "Donor 5",
    "CMO310": "Donor 5",
    "CMO311": "Donor 6",
    "CMO312": "Donor 6"
}

data_dict = {
    "Ch7": {
        "Illumina": "demuxEM/Ch7_scRNA_Illumina_demux.zarr.zip",
        "AVITI": "demuxEM/Ch7_scRNA_AVITI_demux.zarr.zip",
    },
}

def gen_demux_singlet_plot(illu_file, aviti_file, out_prefix):
    data_illu = pg.read_input(illu_file)
    pg.qc_metrics(data_illu, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_illu)

    data_aviti = pg.read_input(aviti_file)
    pg.qc_metrics(data_aviti, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_aviti)

    # Plot regarding CMO tags
    cell_dist_illu = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'assignment'].astype('str').value_counts()
    cell_dist_aviti = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'assignment'].astype('str').value_counts()
    cell_dist_illu = cell_dist_illu.sort_index(ascending=True)
    cell_dist_aviti = cell_dist_aviti.sort_index(ascending=True)
    df_plot = pd.DataFrame({
        'hashtag': cell_dist_illu.index.tolist() + cell_dist_aviti.index.tolist(),
        'counts': cell_dist_illu.values.tolist() + cell_dist_aviti.values.tolist(),
        'technology': ['Illumina'] * cell_dist_illu.size + ['AVITI'] * cell_dist_aviti.size,
    })

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=40)
    sns.barplot(data=df_plot, x='hashtag', y='counts', hue='technology', ax=ax)
    fig.savefig(f"{out_prefix}.demux_singlet.hashtag.png", dpi=500)
    plt.close()

    # Plot regarding Donors
    data_illu.obs['donor'] = ''
    data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'] = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: assignment2donor[s])
    data_aviti.obs['donor'] = ''
    data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'] = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: assignment2donor[s])

    cell_dist_illu = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'].astype('str').value_counts()
    cell_dist_aviti = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'].astype('str').value_counts()
    cell_dist_illu = cell_dist_illu.sort_index(ascending=True)
    cell_dist_aviti = cell_dist_aviti.sort_index(ascending=True)

    df_plot = pd.DataFrame({
        'donor': cell_dist_illu.index.tolist() + cell_dist_aviti.index.tolist(),
        'counts': cell_dist_illu.values.tolist() + cell_dist_aviti.values.tolist(),
        'technology': ['Illumina'] * cell_dist_illu.size + ['AVITI'] * cell_dist_aviti.size,
    })

    fig, ax = plt.subplots(figsize=(10, 8))
    #ax.set_xticklabels(ax.get_xticklabels(), rotation=40)
    sns.barplot(data=df_plot, x='donor', y='counts', hue='technology', ax=ax)
    fig.savefig(f"{out_prefix}.demux_singlet.donor.png", dpi=500)
    plt.close()


if __name__ == "__main__":
    for k, v in data_dict.items():
        print(f"Processing for {k}...")
        gen_demux_singlet_plot(v["Illumina"], v["AVITI"], k)
