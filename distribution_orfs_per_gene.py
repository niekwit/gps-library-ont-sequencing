import matplotlib.pyplot as plt
import pandas as pd

# load original annotations
ori_anno = pd.read_csv("uORF_ORF81.csv")

# rename Gene_ID column to gene_name, and IOH_ID column to orf_id for easier merging
ori_anno = ori_anno.rename(columns={"Gene_ID": "gene_name", "IOH_ID": "orf_id"})

# keep only the relevant columns for merging with the ONT data
ori_anno = ori_anno[["gene_name", "orf_id"]].drop_duplicates()
ori_anno["annotation_source"] = "original"

# Load the consensus ORF data
ont_anno = pd.read_csv(
    "/media/niek/4TB_SSD2/analyses/GPS_ONT/orfs_consensus_with_orf_ids.csv"
)

# from ONT-derived data, keep only the relevant columns for merging with the original annotations
ont_anno = ont_anno[["gene_name", "orf_id"]].drop_duplicates()
ont_anno["annotation_source"] = "ONT-derived"

# combine the two datasets
combined_anno = pd.concat([ori_anno, ont_anno], ignore_index=True)

# count the number of ORFs per gene for each annotation source
orf_counts = (
    combined_anno.groupby(["gene_name", "annotation_source"])
    .size()
    .reset_index(name="orf_count")
)
# plot the distribution of ORFs per gene for each annotation source
# limit x axis to 15
plt.figure(figsize=(12, 6))
plt.xlim(0, 15)
for source in orf_counts["annotation_source"].unique():
    subset = orf_counts[orf_counts["annotation_source"] == source]
    plt.hist(
        subset["orf_count"],
        bins=range(1, subset["orf_count"].max() + 2),
        alpha=0.5,
        label=source,
    )
plt.xlabel("Number of ORFs per Gene")
plt.ylabel("Frequency")
plt.title("Distribution of ORFs per Gene by Annotation Source")
plt.legend()
plt.savefig("orfs_per_gene_distribution.png")
