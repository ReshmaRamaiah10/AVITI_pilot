import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib.patches import Patch


files = {
    "Ch10 BCR": {
        "Illumina": {
            "R1": "quality_score_csvs/Ch10_BCR_ILLUMINA_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_BCR_ILLUMINA_R2_QC.csv",
        },
        "AVITI": {
            "R1": "quality_score_csvs/Ch10_BCR_AVITI_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_BCR_AVITI_R2_QC.csv",
        }
    },
    "Ch10 TCR": {
        "Illumina": {
            "R1": "quality_score_csvs/Ch10_TCR_ILLUMINA_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_TCR_ILLUMINA_R2_QC.csv",
        },
        "AVITI": {
            "R1": "quality_score_csvs/Ch10_TCR_AVITI_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_TCR_AVITI_R2_QC.csv",
        }
    },
    "Ch10": {
        "Illumina": {
            "R1": "quality_score_csvs/Ch10_scRNA_ILLUMINA_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_scRNA_ILLUMINA_R2_QC.csv",
        },
        "AVITI": {
            "R1": "quality_score_csvs/Ch10_scRNA_AVITI_R1_QC.csv",
            "R2": "quality_score_csvs/Ch10_scRNA_AVITI_R2_QC.csv",
        }
    },
    "Ch7": {
        "Illumina": {
            "R1": "quality_score_csvs/Ch7_scRNA_ILLUMINA_R1_QC.csv",
            "R2": "quality_score_csvs/Ch7_scRNA_ILLUMINA_R2_QC.csv",
        },
        "AVITI": {
            "R1": "quality_score_csvs/Ch7_scRNA_AVITI_R1_QC.csv",
            "R2": "quality_score_csvs/Ch7_scRNA_AVITI_R2_QC.csv",
        }
    },
    "Ch7 CellPlex": {
        "Illumina": {
            "R1": "quality_score_csvs/Ch7_CellPlex_ILLUMINA_R1_QC.csv",
            "R2": "quality_score_csvs/Ch7_CellPlex_ILLUMINA_R2_QC.csv",
        },
        "AVITI": {
            "R1": "quality_score_csvs/Ch7_CellPlex_AVITI_R1_QC.csv",
            "R2": "quality_score_csvs/Ch7_CellPlex_AVITI_R2_QC.csv",
        }
    }
}

# Position starts from 0.
r1_regions = {
    "CBC": (0, 15, "#FE8B04"),
    "UMI": (16, 25, "#1B88CC"),
}

r2_regions = {
    "cDNA": (0, 89, "#19384F")
}


def plot_bar(sample_type, filename, ax, region_map, cax, bar_width=1):
    if filename is None:
        return

    df_orig = pd.read_csv(filename, header=None)
    if sample_type in ["Illumina, R1", "AVITI, R1"]:
        df_orig = df_orig[0:26].copy()
    else:
        df_orig = df_orig[0:89].copy()

    bin_count1 = df_orig.loc[:, 40:50].sum(axis=1).values
    bin_count2 = df_orig.loc[:, 30:39].sum(axis=1).values
    bin_count3 = df_orig.loc[:, 20:29].sum(axis=1).values
    bin_count4 = df_orig.loc[:, 10:19].sum(axis=1).values
    bin_count5 = df_orig.loc[:, 0:9].sum(axis=1).values

    df = pd.DataFrame({"40-50": bin_count1, "30-39": bin_count2, "20-29": bin_count3, "10-19": bin_count4, "0-9": bin_count5})

    # Convert to proportion
    s = df.sum(axis=1)
    df = df.apply(lambda col: col / s)

    df.plot(
        kind='bar',
        stacked=True,
        legend=False,
        color=[
            "#818c94",  # 40-50
            "#D3D7DA",  # 30-39
            "#17A8A0",  # 20-29
            "#A3E4AE",  # 10-19
            "#E498CB",  # 0-9
        ],
        width=bar_width,
        ax=ax,
    )
    ax.set_title(sample_type)
    ax.get_xaxis().set_visible(False)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 25, 50, 75, 100])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if sample_type in ["Illumina, R1", "AVITI, R1"]:
        xlimit = 25
        cax_ticks = [-1, 9, 19]
        cax_tick_labels = [0, 10, 20]
    elif sample_type in ["Illumina, R2", "AVITI, R2"]:
        xlimit = 89
        cax_ticks = [-1, 24, 49, 74]
        cax_tick_labels = [0, 25, 50, 75]

    ax.set_xlim(-1-0.5*bar_width, xlimit+0.5*bar_width)


    if len(region_map.keys()) > 1:
        vline_pos = [pos+0.5*bar_width for _, pos, _ in region_map.values()]
        ax.vlines(x=vline_pos[:-1], ymin=0, ymax=1, ls='dashed', color='black')

    # Plot base regions
    cax.get_yaxis().set_visible(False)
    for bar_start, bar_end, color in region_map.values():
        rect = patches.Rectangle(xy=(bar_start-0.5*bar_width, 0), width=(bar_end-bar_start+1)*bar_width, height=1, color=color)
        cax.add_patch(rect)
    if len(region_map.keys()) > 1:
        cax.vlines(x=vline_pos[:-1], ymin=0, ymax=1, ls='dashed', color='black')

    cax.set_xlim(-1-0.5*bar_width, xlimit+0.5*bar_width)
    cax.set_xticks(cax_ticks)
    cax.set_xticklabels(cax_tick_labels, rotation=0)
    cax.spines["left"].set_visible(False)
    cax.spines["top"].set_visible(False)
    cax.spines["right"].set_visible(False)


for sample in files.keys():
    fig, axes = plt.subplots(2, 4, figsize=(20,5), gridspec_kw={'height_ratios': [95, 5], 'hspace': 0.02}, sharey=True, sharex=False)
    plot_bar("Illumina, R1", files[sample]["Illumina"]["R1"], axes[0][0], r1_regions, axes[1][0])
    plot_bar("AVITI, R1", files[sample]["AVITI"]["R1"], axes[0][1], r1_regions, axes[1][1])
    plot_bar("Illumina, R2", files[sample]["Illumina"]["R2"], axes[0][2], r2_regions, axes[1][2])
    plot_bar("AVITI, R2", files[sample]["AVITI"]["R2"], axes[0][3], r2_regions, axes[1][3])


    # Legend for quality bins
    axes[0][-1].legend(loc="right", bbox_to_anchor=(1.5, 0.8), title="Quality bin")

    # Legend for regions
    region_legend_elements = list()
    for region_name, region_setting in (r1_regions|r2_regions).items():
        region_legend_elements.append(Patch(facecolor=region_setting[2], label=region_name))

    axes[1][-1].legend(handles=region_legend_elements, loc="right", bbox_to_anchor=(1.5, 10), title="Region")

    plt.minorticks_off()
    axes[0][0].set_ylabel("Percent reads")
    fig.suptitle(f"{sample} GEX libraries (5')")
    fig.supxlabel("Position in read (bases)")
    fig.savefig(f"figure/{sample}.gex.pair_end_after_trim.png", dpi=500)
    plt.close()
