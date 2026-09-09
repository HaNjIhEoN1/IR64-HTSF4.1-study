#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch, FancyArrowPatch

# ============================================================
# Input files
# ============================================================
BED_PATH = "IR64HTSF4_origin.bed"
FAI_PATH = "IR64-htsf4_pseudo.fasta.fai"

# Output prefix
OUT_PREFIX = "IR64HTSF4_N22_heatQTL"

# N22-derived segment color
SEGMENT_COLOR = "#960018"

# Figure text sizes
QTL_LABEL_FONTSIZE = 9
CHROM_LABEL_FONTSIZE = 11
COORD_LABEL_FONTSIZE = 8
LEGEND_FONTSIZE = 12

# ============================================================
# Heat-stress QTL / marker coordinates on IR64-HTSF4 assembly
# chrom, start, end, label
# ============================================================
QTLS = [
    ("IR64-htsf4_Chr03", 25682591, 26730673, "qSSIY3.1"),
    ("IR64-htsf4_Chr04", 16836792, 16836992, "qHTSF4.1"),
    ("IR64-htsf4_Chr05", 26783016, 27135843, "qSTIY5.1/qSSIY5.1"),
    ("IR64-htsf4_Chr06", 27850539, 27850785, "qSSR6-1"),
    ("IR64-htsf4_Chr07", 29387467, 29387569, "qSSR7-1"),
    ("IR64-htsf4_Chr08", 22914067, 22914208, "qSSR8-1"),
    ("IR64-htsf4_Chr09", 18243760, 18650214, "qSTIPSS9.1"),
    ("IR64-htsf4_Chr09", 20366567, 20366791, "qSSR9-1"),
    ("IR64-htsf4_Chr11", 30505572, 30505728, "qSSR11-1"),
    ("IR64-htsf4_Chr12", 8193928, 9160792, "qSSIPSS12.1"),
]

# Vertical position adjustment.
# qSSR9-1 is intentionally placed between the genome track and qSTIPSS9.1.
TRACK_OFFSET = {
    "qSSR9-1": 0.285,
    "qSTIPSS9.1": 0.445,
}

DEFAULT_TRACK_OFFSET = 0.31

# Very short loci are expanded visually only.
# Their actual genomic coordinate is NOT changed.
SHORT_LOCUS_THRESHOLD = 100_000
SHORT_LOCUS_VISIBLE_FRACTION = 0.015


def read_inputs():
    bed = pd.read_csv(
        BED_PATH,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "origin"]
    )

    fai = pd.read_csv(
        FAI_PATH,
        sep="\t",
        header=None,
        names=["chrom", "length", "offset", "linebases", "linewidth"]
    )

    return bed, fai


def validate_inputs(bed, fai):
    """
    Three-stage input validation.
    """

    # --------------------------------------------------------
    # Validation 1: file/schema
    # --------------------------------------------------------
    required_bed = {"chrom", "start", "end", "origin"}
    if not required_bed.issubset(bed.columns):
        raise ValueError("BED schema is invalid.")

    if len(fai) != 12:
        raise ValueError(f"Expected 12 chromosomes in FAI, found {len(fai)}.")

    origins = set(bed["origin"].dropna().unique())
    if not origins.issubset({"IR64", "N22"}):
        raise ValueError(f"Unexpected origin labels: {origins}")

    # --------------------------------------------------------
    # Validation 2: BED coordinates
    # --------------------------------------------------------
    lengths = dict(zip(fai["chrom"], fai["length"]))

    if not bed["chrom"].isin(lengths).all():
        bad = sorted(set(bed.loc[~bed["chrom"].isin(lengths), "chrom"]))
        raise ValueError(f"BED chromosomes absent from FAI: {bad}")

    if (bed["start"] < 0).any():
        raise ValueError("Negative BED start coordinate detected.")

    if (bed["end"] <= bed["start"]).any():
        raise ValueError("BED row with end <= start detected.")

    for row in bed.itertuples():
        if row.end > lengths[row.chrom]:
            raise ValueError(
                f"BED coordinate exceeds chromosome length: "
                f"{row.chrom}:{row.start}-{row.end}"
            )

    # --------------------------------------------------------
    # Validation 3: QTL coordinates
    # --------------------------------------------------------
    for chrom, start, end, name in QTLS:
        if chrom not in lengths:
            raise ValueError(f"{name}: chromosome {chrom} absent from FAI.")

        if not (0 <= start < end <= lengths[chrom]):
            raise ValueError(
                f"{name}: invalid coordinate {chrom}:{start}-{end}"
            )

    return lengths


def merge_n22_blocks(bed):
    """
    Merge overlapping or directly adjacent N22 windows.
    """
    n22 = (
        bed.loc[bed["origin"] == "N22"]
        .sort_values(["chrom", "start", "end"])
        .copy()
    )

    merged = []

    for chrom, group in n22.groupby("chrom", sort=False):
        current_start = None
        current_end = None

        for row in group.itertuples():
            if current_start is None:
                current_start = row.start
                current_end = row.end

            elif row.start <= current_end:
                current_end = max(current_end, row.end)

            else:
                merged.append((chrom, current_start, current_end))
                current_start = row.start
                current_end = row.end

        if current_start is not None:
            merged.append((chrom, current_start, current_end))

    return pd.DataFrame(
        merged,
        columns=["chrom", "start", "end"]
    )


def short_chr_name(chrom):
    return chrom.replace("IR64-htsf4_", "")


def draw_qtl_track(ax, chr_len, start, end, name, chromosome_y):
    """
    Draw QTL/marker above chromosome.

    The arrow points DOWNWARD toward the genome track.

    Long QTL:
        actual genomic interval is drawn as an open box.

    Short marker-level locus:
        genomic midpoint is retained, but the box width is enlarged
        schematically so that it remains visible in a genome-wide figure.
    """

    chr_height = 0.24
    box_height = 0.085

    genome_top_y = chromosome_y + chr_height / 2
    track_y = chromosome_y + TRACK_OFFSET.get(
        name,
        DEFAULT_TRACK_OFFSET
    )

    locus_width = end - start

    if locus_width < SHORT_LOCUS_THRESHOLD:
        midpoint = (start + end) / 2

        visible_width = chr_len * SHORT_LOCUS_VISIBLE_FRACTION
        x0 = max(0, midpoint - visible_width / 2)
        x1 = min(chr_len, midpoint + visible_width / 2)

    else:
        x0 = start
        x1 = end
        midpoint = (x0 + x1) / 2

    # Open QTL/marker box
    ax.add_patch(
        Rectangle(
            (x0, track_y - box_height / 2),
            x1 - x0,
            box_height,
            fill=False,
            linewidth=1.0
        )
    )

    # Genome-facing arrow
    ax.add_patch(
        FancyArrowPatch(
            (midpoint, track_y - box_height / 2),
            (midpoint, genome_top_y + 0.01),
            arrowstyle="-|>",
            mutation_scale=7.5,
            linewidth=1.0
        )
    )

    return midpoint, track_y


def plot_genome(
    bed,
    fai,
    lengths,
    merged_n22,
    output_file,
    show_labels=True
):

    chrom_order = fai["chrom"].tolist()
    max_len = fai["length"].max()

    fig, ax = plt.subplots(figsize=(13.5, 9.0))

    y_positions = {
        chrom: len(chrom_order) - 1 - i
        for i, chrom in enumerate(chrom_order)
    }

    chr_height = 0.24

    for chrom in chrom_order:

        y = y_positions[chrom]
        chr_len = lengths[chrom]

        # ----------------------------------------------------
        # Chromosome backbone
        # ----------------------------------------------------
        ax.add_patch(
            Rectangle(
                (0, y - chr_height / 2),
                chr_len,
                chr_height,
                fill=False,
                linewidth=0.9
            )
        )

        # ----------------------------------------------------
        # N22-derived segments
        # ----------------------------------------------------
        chr_n22 = merged_n22.loc[
            merged_n22["chrom"] == chrom
        ]

        for row in chr_n22.itertuples():
            ax.broken_barh(
                [(row.start, row.end - row.start)],
                (y - chr_height / 2, chr_height),
                facecolors=SEGMENT_COLOR
            )

        # ----------------------------------------------------
        # Heat-stress QTL / markers
        # ----------------------------------------------------
        for qchrom, start, end, name in QTLS:

            if qchrom != chrom:
                continue

            label_x, label_y = draw_qtl_track(
                ax,
                chr_len,
                start,
                end,
                name,
                y
            )

            if show_labels:
                ax.annotate(
                    name,
                    (label_x, label_y + 0.05),
                    xytext=(0, 2),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=QTL_LABEL_FONTSIZE
                )

        # ----------------------------------------------------
        # Chromosome label
        # ----------------------------------------------------
        ax.text(
            -0.025 * max_len,
            y,
            short_chr_name(chrom),
            ha="right",
            va="center",
            fontsize=CHROM_LABEL_FONTSIZE
        )

        # ----------------------------------------------------
        # 10-Mb coordinate ticks
        # ----------------------------------------------------
        tick = 0

        while tick <= chr_len:

            ax.vlines(
                tick,
                y - chr_height / 2 - 0.035,
                y - chr_height / 2,
                linewidth=0.6
            )

            ax.text(
                tick,
                y - chr_height / 2 - 0.08,
                f"{tick / 1e6:.0f}",
                ha="center",
                va="top",
                fontsize=COORD_LABEL_FONTSIZE
            )

            tick += 10_000_000

    # --------------------------------------------------------
    # Axes
    # --------------------------------------------------------
    ax.set_xlim(
        -0.06 * max_len,
        max_len * 1.02
    )

    ax.set_ylim(
        -0.7,
        len(chrom_order) + 0.35
    )

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    ax.legend(
        handles=[
            Patch(
                facecolor=SEGMENT_COLOR,
                edgecolor=SEGMENT_COLOR,
                label="N22-derived segment"
            ),
            Patch(
                fill=False,
                edgecolor="black",
                label="Heat-stress QTL / marker"
            )
        ],
        loc="lower right",
        frameon=False,
        fontsize=LEGEND_FONTSIZE
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=600,
        bbox_inches="tight"
    )

    plt.close(fig)


def main():

    bed, fai = read_inputs()

    lengths = validate_inputs(
        bed,
        fai
    )

    merged_n22 = merge_n22_blocks(
        bed
    )

    print(
        f"[INFO] BED rows: {len(bed):,}"
    )

    print(
        f"[INFO] N22 windows: "
        f"{(bed['origin'] == 'N22').sum():,}"
    )

    print(
        f"[INFO] Merged N22 blocks: "
        f"{len(merged_n22):,}"
    )

    print(
        f"[INFO] QTL/marker loci: "
        f"{len(QTLS)}"
    )

    # clean
    plot_genome(
        bed,
        fai,
        lengths,
        merged_n22,
        f"{OUT_PREFIX}_clean.png",
        show_labels=False
    )

    # labeled PNG
    plot_genome(
        bed,
        fai,
        lengths,
        merged_n22,
        f"{OUT_PREFIX}_labeled.png",
        show_labels=True
    )

    # labeled vector PDF
    plot_genome(
        bed,
        fai,
        lengths,
        merged_n22,
        f"{OUT_PREFIX}_labeled.pdf",
        show_labels=True
    )

    print("[DONE]")


if __name__ == "__main__":
    main()


