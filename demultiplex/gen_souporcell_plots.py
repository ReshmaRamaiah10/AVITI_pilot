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


def gen_ch7_plots():
    # Load Souporcell results
    data_illu = pg.read_input("souporcell/Ch7_scRNA_Illumina_demux.zarr.zip")
    pg.qc_metrics(data_illu, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_illu)

    data_aviti = pg.read_input("souporcell/Ch7_scRNA_AVITI_demux.zarr.zip")
    pg.qc_metrics(data_aviti, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_aviti)

    # Rename Souporcell assignment to avoid confusion
    soup_remap = {
        'Donor1': 'CB1',
        'Donor2': 'CB2',
        'Donor3': 'CB3',
        'Donor4': 'CB4',
        'Donor5': 'CB5',
        'Donor6': 'CB6',
    }

    data_aviti.obs['donor'] = ''
    data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'] = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: soup_remap[s])

    data_illu.obs['donor'] = ''
    data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'] = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: soup_remap[s])

    # Load DemuxEM results as reference
    data_illu_em = pg.read_input("demuxEM/Ch7_scRNA_Illumina_demux.zarr.zip")
    pg.qc_metrics(data_illu_em, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_illu_em)

    data_aviti_em = pg.read_input("demuxEM/Ch7_scRNA_AVITI_demux.zarr.zip")
    pg.qc_metrics(data_aviti_em, min_genes=500, mito_prefix="MT-", percent_mito=20)
    pg.filter_data(data_aviti_em)

    cmo2donor = {
        'CMO301': 'Donor 1',
        'CMO302': 'Donor 1',
        'CMO303': 'Donor 2',
        'CMO304': 'Donor 2',
        'CMO305': 'Donor 3',
        'CMO306': 'Donor 3',
        'CMO307': 'Donor 4',
        'CMO308': 'Donor 4',
        'CMO309': 'Donor 5',
        'CMO310': 'Donor 5',
        'CMO311': 'Donor 6',
        'CMO312': 'Donor 6',
    }

    data_aviti_em.obs['donor'] = ''
    data_aviti_em.obs.loc[data_aviti_em.obs['demux_type']=='singlet', 'donor'] = data_aviti_em.obs.loc[data_aviti_em.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: cmo2donor[s])

    data_illu_em.obs['donor'] = ''
    data_illu_em.obs.loc[data_illu_em.obs['demux_type']=='singlet', 'donor'] = data_illu_em.obs.loc[data_illu_em.obs['demux_type']=='singlet', 'assignment'].astype(str).map(lambda s: cmo2donor[s])

    idx_singlet_aviti = np.intersect1d(data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet'].index, data_aviti_em.obs.loc[data_aviti_em.obs['demux_type']=='singlet'].index)
    idx_singlet_illu = np.intersect1d(data_illu.obs.loc[data_illu.obs['demux_type']=='singlet'].index, data_illu_em.obs.loc[data_illu_em.obs['demux_type']=='singlet'].index)

    df_aviti = pd.crosstab(data_aviti_em.obs.loc[idx_singlet_aviti, 'donor'], data_aviti.obs.loc[idx_singlet_aviti, 'donor'])
    df_illu = pd.crosstab(data_illu_em.obs.loc[idx_singlet_illu, 'donor'], data_illu.obs.loc[idx_singlet_illu, 'donor'])

    # Plot confusion matrix
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(data=df_aviti, annot=True)
    fig.savefig("ch7_aviti_confusion.png", dpi=500)
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(data=df_illu, annot=True)
    fig.savefig("ch7_illu_confusion.png", dpi=500)
    plt.close()

    # Remap Souporcell assignment due to confusion matrix
    aviti_donor_map = {
        'CB1': 'Donor 2',
        'CB2': 'Donor 4',
        'CB3': 'Donor 3',
        'CB4': 'Donor 1',
        'CB5': 'Donor 5',
        'CB6': 'Donor 6',
    }
    data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'] = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'].map(lambda s: aviti_donor_map[s])

    illu_donor_map = {
        'CB1': 'Donor 2',
        'CB2': 'Donor 1',
        'CB3': 'Donor 6',
        'CB4': 'Donor 3',
        'CB5': 'Donor 4',
        'CB6': 'Donor 5',
    }
    data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'] = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'].map(lambda s: illu_donor_map[s])

    # Plot singlet plot
    cell_dist_illu = data_illu.obs.loc[data_illu.obs['demux_type']=='singlet', 'donor'].astype(str).value_counts()
    cell_dist_aviti = data_aviti.obs.loc[data_aviti.obs['demux_type']=='singlet', 'donor'].astype(str).value_counts()
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
    fig.savefig("Ch7.souporcell_singlet.png", dpi=500)
    plt.close()


if __name__ == "__main__":
    gen_ch7_plots()
