import sys
import pysam
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import Counter


def get_barcode_counts(bam_path):
    """Extracts BC tags and returns a sorted DataFrame of counts."""
    counts = Counter()

    print(f"Reading BAM: {bam_path}...")
    total_reads = pysam.AlignmentFile(bam_path, "rb").mapped
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in tqdm(bam, desc="Processing reads", total=total_reads):
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Check if the read has the BC tag
            if read.has_tag("BC"):
                barcode = read.get_tag("BC")
                counts[barcode] += 1

    # Convert to DataFrame
    df = pd.DataFrame.from_dict(counts, orient="index", columns=["count"])
    df.index.name = "barcode"
    return df.reset_index()


def plot_bam_knee(bam_path, unique):
    # Get counts from BAM
    df = get_barcode_counts(bam_path)

    # Data Prep: Sort and Rank
    df = df.sort_values(by="count", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1

    sns.set_context("notebook", font_scale=1.1)
    sns.set_style("ticks")

    fig, ax = plt.subplots(figsize=(10, 6))

    # Main Knee Plot
    sns.lineplot(data=df, x="rank", y="count", color="navy", linewidth=1.5, ax=ax)

    # Log-Log Scale
    ax.set(xscale="log", yscale="log")

    # Set tight limits to prevent "overhang" or empty space
    # round unique up to nearest 10K for better visualization
    unique_limit = ((unique + 9999) // 10000) * 10000
    ax.set_xlim(df["rank"].min(), unique_limit)
    ax.set_ylim(0.8, df["count"].max() + 10)

    # Reference Lines
    ax.axvline(
        x=180000, color="crimson", linestyle="--", label="Expected (180k)", alpha=0.8
    )
    # Convert unique to XXXk format for label
    unique_k = unique / 1000
    ax.axvline(
        x=unique,
        color="forestgreen",
        linestyle=":",
        label=f"Current ({unique_k:.0f}k)",
        alpha=0.8,
    )

    # offset=10 creates the gap at (0,0)
    # trim=False ensures the lines are solid/continuous at the ends
    sns.despine(ax=ax, offset=10, trim=False)

    # Labels and Legend
    ax.set_xlabel("Barcode Rank (Log10)")
    ax.set_ylabel("Reads per Barcode (Log10)")
    ax.legend(frameon=False)

    # Add title with the BAM filename (without path and extension)
    bam_name = bam_path.split("/")[-1].replace(".bam", "")
    ax.set_title(f"Knee Plot for {bam_name}", fontsize=14)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Get command line arguments
    args = sys.argv[1:]
    if len(args) != 2:
        print("Usage: python script.py <bam_file> <limit for x axis>")
        sys.exit(1)

    bam = args[0]
    unique = int(args[1])

    plot_bam_knee(bam, unique)
