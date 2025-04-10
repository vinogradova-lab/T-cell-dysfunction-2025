import textwrap
import json
import requests
import math
import re
import pandas as pd
import numpy as np


import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import logomaker as logomaker

# constants for domain plot
LINE_WIDTH = 0.5
PLT_DIMENSIONS = (10, 5)

# pagination size for uniprot queries
CHUNK_SIZE = 500


def create_entry_cache(df):
    """Given a list of cysteine aggregation files,
    create a cache of UniProt entries with each unique
    protein.
    """
    entry_dict = dict()
    uniprots = set()  # use set so each identifier is unique
    for uniprot in df["uniprot"]:
        if "contaminant" not in uniprot:
            if "reversed" not in uniprot:
                uniprots.add(uniprot)

    # Break the list of ids into smaller lists to not overwhelm uniprot
    i = 0
    uniprots = list(uniprots)  # convert to list so we can subscript
    chunks = [uniprots[x : x + CHUNK_SIZE] for x in range(0, len(uniprots), CHUNK_SIZE)]
    print("Querying UniProt...")
    for chunk in chunks:
        chunk_num = i * CHUNK_SIZE
        print("Retrieved " + str(chunk_num) + " entries out of " + str(len(uniprots)))
        i = i + 1
        entries = UniprotkbClient.fetch_many(chunk)
        for entry in entries:
            entry_dict[entry["primaryAccession"]] = entry
    print("Done querying UniProt.")
    return entry_dict


def get_features(uniprot, cys, entry_dict):
    """Given a uniprot identifier, a cysteine residue,
    and a cache of UniProt entries, retrieve all features
    matched to the cysteine residue.
    """
    ret_list = list()
    if uniprot not in entry_dict:
        return ret_list
    entry = entry_dict[uniprot]
    cys_position = int(cys.split("C")[1])
    if "features" not in entry:
        return ret_list
    for feature in entry["features"]:
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
        if (
            type(start) == int
            and type(end) == int
            and cys_position >= start
            and cys_position <= end
        ):
            # Nearly every protein has a chain annotation for the entire length of the sequence
            # these are not informative so drop them.

            if (not feature["type"] == "Chain") and (feature["type"] == "Domain"):
                ret_list.append(feature)

    return ret_list


def get_feature_rows(residue, row, entry_dict):
    """Given a residue and a cache of UniProt entries,
    return a list of rows corresponding to feature-residue
    matches. If the residue does not match to any features,
    return a row with None labels.
    """
    fs = get_features(row["uniprot"], residue, entry_dict)
    data = []
    none = "None"
    if len(fs) == 0:
        row["feature_description"] = none
        row["feature_type"] = none
        data.append(row)
    else:
        for feature in fs:
            # Parse the feature dict into variables for our row
            feature_row = row.copy()
            feature_row["feature_description"] = feature["description"]
            feature_row["feature_type"] = feature["type"]
            feature_row["feature_start"] = feature["location"]["start"]["value"]
            feature_row["feature_end"] = feature["location"]["end"]["value"]
            evidences = list()
            evidence_codes = list()
            # List the sources of evidence for this feature
            if "evidences" in feature:
                for evidence in feature["evidences"]:
                    # Some features do not report their source
                    if "source" in evidence:
                        evidences.append(evidence["source"])
                    if "evidenceCode" in evidence:
                        evidence_codes.append(evidence["evidenceCode"])
                feature_row["feature_evidence_source"] = ", ".join(evidences)
                feature_row["feature_evidence_code"] = ", ".join(evidence_codes)
            data.append(feature_row)
    return data


def create_feature_df(df, entry_dict):
    """Given a path to a Cysteine aggregation output file,
    and optionally a cache of UniProt entries,
    annotate each residue with functional domain information
    """
    print("Matching residues to features..")
    data = []
    for i, row in df.iterrows():
        # Some rows have multiple residues. Create a separate row
        # for each residue
        residues = re.split(";|,", row["residue"])  # .split(r",|;")
        row = row.copy()
        for residue in residues:
            row["residue"] = residue
            rows = get_feature_rows(residue, row, entry_dict)
            for row in rows:
                data.append(row)
    return pd.DataFrame(data)


def get_protein_length(uniprot):
    url = "https://rest.uniprot.org/uniprotkb/{}?format=json".format(uniprot)
    uniprot_domains = requests.get(url).text
    return json.loads(uniprot_domains)["sequence"]["length"]


def plot_domains(uniprot, ax, protein_range, residues, domains_dict = {}):
    if len(domains_dict) == 0:
        # retrieve domain info from uniprot
        url = "https://rest.uniprot.org/uniprotkb/{}?format=json&fields=ft_domain,cc_domain,ft_region,ft_motif".format(
            uniprot
        )
        uniprot_domains = requests.get(url).text
        domains_dict = {}
        for domain in json.loads(uniprot_domains)["features"]:
            name = domain["description"]
            start = domain["location"]["start"]["value"]
            end = domain["location"]["end"]["value"]
            domains_dict[name] = [start, end]

    # plot domains

    with sns.plotting_context("paper", font_scale=0.6):
        height = 1.0
        # background patch
        ax.add_patch(
            patches.FancyBboxPatch(
                (protein_range[0], 0.5 - height / 2),
                protein_range[1] - protein_range[0],
                height,
                boxstyle=f"round4, pad=0.0, rounding_size=0.0",
                linewidth=LINE_WIDTH,
                edgecolor="black",
                facecolor="white",
                alpha=1.0,
                zorder=1,
            )
        )
        for name, (start, end) in domains_dict.items():
            length = end - start
            ax.add_patch(
                patches.FancyBboxPatch(
                    (start, 0.5 - height / 2),
                    length,
                    height,
                    boxstyle=f"round4, pad=0.0, rounding_size=0.0",
                    linewidth=LINE_WIDTH,
                    edgecolor="black",
                    facecolor="#D3D3D3",
                    alpha=1.0,
                    zorder=1,
                )
            )
            domain_name = textwrap.fill(name, 10)
            ax.text(
                start + length / 2,
                0.5,
                domain_name,
                ha="center",
                va="center",
                fontsize=10,
            )

        # ax.axhline(y=0.5, xmin=0, xmax=5000, linewidth=1, color="black", zorder=0)
        ax.spines[["right", "top", "left", "bottom"]].set_visible(False)
        ax.set_xlim(protein_range[0], protein_range[1])
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xticks(residues["pos"])
        ax.set_xticklabels(residues["residue"], fontdict={"weight": "bold"})
        tick_labels = ax.get_xticklabels()
        for i, hue in enumerate(residues["hue"]):
            tick_labels[i].set_color(hue)


def plot_reactivity_changes(uniprot, condition, axs, protein_range, df, scatter_data):
    df = df[df["condition"] == condition]

    plot_palette = {"#404040": "#404040", "blue": "blue", "red": "red"}
    plot_alpha = 0.5
    sns.boxplot(
        data=df,
        x="pos",
        y="LFC_tmt_abpp",
        hue="hue",
        palette=plot_palette,
        order=range(protein_range[0], protein_range[1]),
        boxprops=dict(alpha=plot_alpha),
        showfliers=False,
        linewidth=LINE_WIDTH,
        ax=axs,
        width=15,
    )
    sns.swarmplot(
        data=df,
        x="pos",
        y="LFC_tmt_abpp",
        hue="hue",
        palette=plot_palette,
        linewidth=LINE_WIDTH,
        alpha=1,
        edgecolor="black",
        order=range(protein_range[0], protein_range[1]),
        ax=axs,
        size=5,
    )
    axs.axhline(y=1, color="black", linestyle="dotted", zorder=0, linewidth=0.5)
    axs.axhline(y=-1, color="black", linestyle="dotted", zorder=0, linewidth=0.5)
    axs.set_ylabel("{}\nLFC".format(condition), fontsize=12)
    axs.legend([], [], frameon=False)
    if condition != "D8C":
        axs.set_xticklabels([])
    limit_data = scatter_data[scatter_data.uniprot == uniprot]
    y_range = [
        min(limit_data["LFC_tmt_abpp"] - 0.1),
        max(limit_data["LFC_tmt_abpp"] + 0.1),
    ]
    axs.set_ylim(y_range)
    axs.set_xlabel("")
    axs.set_xticks(axs.get_xticks()[::50])
    axs.set_yticks(range(math.floor(y_range[0]), math.floor(y_range[1]), 1))
    axs.spines[["right", "top", "left", "bottom"]].set_linewidth(0.5)
    axs.yaxis.set_tick_params(width=LINE_WIDTH, length=2.5)
    axs.xaxis.set_tick_params(width=LINE_WIDTH, length=2.5)


def loc_plot(plot_uniprot, plot_gene, results_dir, scatter_data):

    df = scatter_data.reset_index()
    df = df.loc[df.uniprot == plot_uniprot]
    df.loc[(df["LFC_tmt_abpp_mean"] >= 1), "hue"] = "blue"
    df.loc[(df["LFC_tmt_abpp_mean"] <= -1), "hue"] = "red"
    df.loc[(df["LFC_tmt_abpp_mean"] > -1) & (df["LFC_tmt_abpp_mean"] < 1), "hue"] = (
        "#404040"
    )

    # position of cysteine for plotting
    df = df.assign(
        pos=df.residue.str.split(";", expand=True)[0]
        .str.split("C", expand=True)[1]
        .astype(int)
    ).sort_values("pos", ascending=True)

    protein_range = [0, get_protein_length(plot_uniprot)]
    plt_size = PLT_DIMENSIONS
    fig = plt.figure(figsize=plt_size, num=1, clear=True, dpi=1000)
    # Define the outer 6x1 grid
    gs = matplotlib.gridspec.GridSpec(
        6, 1, figure=fig, height_ratios=[1, 1, 1, 1, 0.25, 0.5], hspace=0.15
    )
    for i, condition in enumerate(["D4A", "D4C", "D8A", "D8C"]):
        plot_reactivity_changes(
            plot_uniprot,
            condition,
            fig.add_subplot(gs[i]),
            protein_range,
            df,
            scatter_data,
        )
    plot_domains(
        uniprot=plot_uniprot,
        ax=fig.add_subplot(gs[-1]),
        protein_range=protein_range,
        residues=df[df["condition"] == "D8C"][
            ["residue", "pos", "hue"]
        ].drop_duplicates(),
    )
    plt.suptitle(plot_gene + " cysteine reactivity")

    plt.show()


def read_fasta(fn):
    seqs = {}
    for record in SeqIO.parse(fn, format="fasta"):
        uniprot_id = str(record.id).split("|")[1]
        seqs[uniprot_id] = str(record.seq).replace("*", "")
    return seqs


def fetch_flanking_seq(r):
    flank_len = 5
    protein_sequence = r["uniprot_sequence"]
    protein_len = len(r["uniprot_sequence"])

    cys_index = r["residue_loc"] - 1

    start = cys_index - flank_len
    if start < 0:
        start = 0

    stop = cys_index + flank_len
    if protein_len < stop:
        stop = protein_len

    r["flanking_seq"] = protein_sequence[start:stop]
    r["flanking_seq_len"] = len(protein_sequence[start:stop])
    return r


def position_frequency_table(seq_df):
    rs = []
    for seq in seq_df["flanking_seq"].str.split(""):
        row = {}
        for i, residue in enumerate(seq):
            row[i - 1] = residue
        rs.append(row)

    p_df = pd.DataFrame(rs).drop([-1, 10], axis=1).melt()
    p_df["count"] = 1
    p_df = p_df.groupby(["variable", "value"]).count().reset_index()

    ls = []
    for pos in set(p_df["variable"]):
        pos_df = p_df[p_df["variable"] == pos].copy()
        pos_df["proportion"] = pos_df["count"] / sum(pos_df["count"])
        ls.append(pos_df)
    return (
        pd.concat(ls)
        .drop("count", axis=1)
        .pivot(columns="value", values="proportion", index="variable")
        .replace(np.nan, 0)
    )


def logo_plot(position_freq_table):
    # create Logo object
    crp_logo = logomaker.Logo(
        position_freq_table,
        shade_below=0.5,
        fade_below=0.5,
        font_name="Arial Rounded MT Bold",
    )

    # style using Logo methods
    crp_logo.style_spines(visible=False)
    crp_logo.style_spines(spines=["left", "bottom"], visible=True)
    crp_logo.style_xticks(fmt="%d", anchor=0, family="Arial", size=6)

    # style using Axes methods
    crp_logo.ax.set_ylabel("frequency", labelpad=-1, family="Arial")
    crp_logo.ax.xaxis.set_ticks_position("none")
    crp_logo.ax.xaxis.set_tick_params(pad=-1)

    # style and show figure
    return crp_logo


def make_logo_plot(ratio_combfiles):
    uniprot = read_fasta(
        "/Users/henrysanford/Dropbox @RU Dropbox/Vinogradova Laboratory/Vinogradova Laboratory/Henry_data processing/03_fasta_files/homo_sapiens/20240312/uniprotkb_proteome_UP000005640_AND_revi_2024_03_12.fasta"
    )
    ratio_combfiles["uniprot_sequence"] = ratio_combfiles["uniprot"].apply(
        lambda x: uniprot[x]
    )
    ratio_combfiles["residue_loc"] = ratio_combfiles["residue"].apply(
        lambda x: int(x.split("C")[1].split(",")[0])
    )
    ratio_combfiles = ratio_combfiles.apply(fetch_flanking_seq, axis=1)
    ratio_combfiles = ratio_combfiles[ratio_combfiles["flanking_seq_len"] == 10]
    pos_freq_table = position_frequency_table(ratio_combfiles)
    return logo_plot(pos_freq_table)
