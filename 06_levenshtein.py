"""
Correct barcode sequences in a BAM file using Levenshtein distance.

1. Setup and Argument Parsing

    Initialization: The script defines parameters for what a "good" barcode looks like. Specifically, it targets sequences of length 24 or 30 (--target-lengths).

    Thresholds: It sets limits for the "correction" (how different a sequence can be from the truth) using Levenshtein distance and a frequency ratio (the "true" barcode must be much more common than the "error").

2. The Main Loop (Grouping by Reference)

    Sequential Processing: The main() function iterates through the BAM file.

    Coordinate-Based Chunking: It identifies when the reference_name changes (e.g., moving from one cDNA transcript to another) and gathers all reads for that specific transcript into a list (reads_in_group).

    Memory Management: By processing one transcript group at a time, it avoids loading the entire file into memory while still allowing for localized barcode comparison.

3. Analyzing Barcodes within a Group (process_group)

    Frequency Counting: It uses Counter to determine how many times each unique barcode sequence appears within that specific genomic location.

    Identifying Centroids: It identifies "Centroid" sequences—these are barcodes that match the user-defined target lengths (24 or 30). These are considered the "potential truths."

    Sorting by Abundance: Unique barcodes are sorted from most frequent to least frequent to prioritize the most likely correct sequences as anchors.

4. The Correction Logic (The "Correction Map")

For every barcode that isn't a perfect centroid, the script searches for a match among the centroids based on four strict criteria:

    Abundance Ratio: The candidate centroid must have at least ratio_threshold times more reads than the current sequence (ensuring we don't collapse two equally frequent sequences).

    Length Window: The length difference between the sequence and the centroid must be within the --window (e.g., ±2 bp).

    Edit Distance: The Levenshtein distance (number of insertions, deletions, or substitutions) must be within --max-edit-dist.

    Best Match Selection: If multiple centroids fit, it chooses the one with the smallest distance; if distances are equal, it chooses the one with the highest read count.

5. Applying Changes and Writing

    Tag Updating: If a correction is found, read.set_tag("BC", new_bc) overwrites the original noisy sequence in the BAM record with the corrected centroid sequence.

    Preservation: Reads without any BC tag or those that don't meet the correction criteria are written to the output file unchanged to preserve data integrity.

6. Metrics and Reporting

    Tracking Stats: The STATS dictionary tracks "Before vs. After" metrics, such as how many unique barcodes were collapsed and the total percentage of reads that were modified.

    Summary: After the file is processed, it prints a report showing the efficiency of the correction—ideally, you want to see the "Unique barcodes (After)" count be significantly lower than the "Before" count.

"""

import argparse
import pysam
from collections import Counter
from rapidfuzz.distance import Levenshtein
import tqdm

# Create command line parser
parser = argparse.ArgumentParser(
    description="Correct barcode sequences in a BAM file using Levenshtein distance."
)
parser.add_argument("-i", "--input", required=True, help="Input BAM file")
parser.add_argument("-o", "--output", required=True, help="Output BAM file")
parser.add_argument(
    "--target-lengths",
    nargs="+",
    type=int,
    default=[24, 30],
    help="Target barcode lengths to consider as centroids (default: 24 30)",
)
parser.add_argument(
    "--window",
    type=int,
    default=2,
    help="Maximum length difference to consider for correction (default: 2)",
)
parser.add_argument(
    "--max-edit-dist",
    type=int,
    default=2,
    help="Maximum Levenshtein distance for correction (default: 2)",
)
parser.add_argument(
    "--ratio-threshold",
    type=int,
    default=5,
    help="Minimum frequency ratio for correction (default: 5)",
)

args = parser.parse_args()

# Global metrics
STATS = {
    "total_reads_with_bc": 0,
    "total_corrected_reads": 0,
    "unique_bcs_before": set(),
    "unique_bcs_after": set(),
}


def process_group(reads, writer):
    # Filter reads that actually have a BC tag
    reads_with_bc = [r for r in reads if r.has_tag("BC")]

    # Write reads without BC immediately and move on
    for r in reads:
        if not r.has_tag("BC"):
            writer.write(r)

    if not reads_with_bc:
        return

    raw_bcs = [r.get_tag("BC") for r in reads_with_bc]

    # Count occurrences of each unique BC tag in the current group
    counts = Counter(raw_bcs)
    unique_bcs = sorted(counts.keys(), key=lambda x: counts[x], reverse=True)

    # Update "Before" stats
    STATS["total_reads_with_bc"] += len(raw_bcs)
    STATS["unique_bcs_before"].update(unique_bcs)

    # Identify "Centroids" (potential true BCs) based on target lengths (24, 30)
    centroids = [b for b in unique_bcs if len(b) in args.target_lengths]
    correction_map = {}

    # For each barcode, find a nearby highly-abundant Centroid to collapse into
    for seq in unique_bcs:
        # If it's already a centroid, it maps to itself and so we can
        # skip the distance calculations
        if seq in centroids:
            correction_map[seq] = seq
            continue

        seq_count = counts[seq]
        candidates = []

        # Compare each non-centroid sequence to all centroids and
        # apply the correction criteria
        for c in centroids:
            # The sequence must have at least ratio_threshold times more
            # reads than the current sequence
            if counts[c] >= (seq_count * args.ratio_threshold):
                # The length difference must be within the specified window
                # This is due to ONT errors often causing small indels, so we
                # allow for some flexibility in length
                if abs(len(seq) - len(c)) <= args.window:
                    # This calculates the minimum number of edits (substitutions,
                    # insertions, deletions) needed to turn the error into the centroid.
                    # It stops calculating if it realizes the distance will definitely
                    # exceed the set limit (e.g. 2).
                    dist = Levenshtein.distance(seq, c, score_cutoff=args.max_edit_dist)

                    # If the distance is within the cutoff, add it to the candidate list
                    if dist <= args.max_edit_dist:
                        candidates.append((dist, counts[c], c))

        # If we found any candidates, choose the one with the
        # smallest distance (and highest count if there's a tie)
        if candidates:
            candidates.sort(key=lambda x: (x[0], -x[1]))
            correction_map[seq] = candidates[0][2]
        else:
            # If no centroid meets the ratio, length, and distance requirements, it
            # is assumed the sequence is unique
            correction_map[seq] = seq

    # Write updated reads and update "After" stats
    for read in reads_with_bc:
        old_bc = read.get_tag("BC")
        new_bc = correction_map.get(old_bc, old_bc)

        if old_bc != new_bc:
            STATS["total_corrected_reads"] += 1

            # Overwrite the BC tag in the BAM record with the corrected sequence
            read.set_tag("BC", new_bc)

        STATS["unique_bcs_after"].add(new_bc)
        writer.write(read)


def main():
    pysam.set_verbosity(0)

    with pysam.AlignmentFile(args.input, "rb") as reader:
        with pysam.AlignmentFile(args.output, "wb", template=reader) as writer:
            current_rname = None
            reads_in_group = []

            for read in tqdm.tqdm(reader, desc="Processing reads"):
                rname = read.reference_name

                # Group reads by reference (cDNA) to process transcript-specific barcodes
                if rname != current_rname:
                    if reads_in_group:
                        process_group(reads_in_group, writer)
                    current_rname = rname
                    reads_in_group = [read]
                else:
                    reads_in_group.append(read)

            if reads_in_group:
                process_group(reads_in_group, writer)

    # Final Summary Report
    print("\n" + "=" * 40)
    print("BARCODE CORRECTION SUMMARY")
    print("=" * 40)
    print(f"Input file:      {args.input}")
    print(f"Target lengths:  {args.target_lengths}")
    print(f"Window size:     {args.window}")
    print(f"Max Edit Dist:   {args.max_edit_dist}")
    print(f"Ratio Threshold: {args.ratio_threshold}")
    print("-" * 40)
    print(f"Total reads processed (with BC):  {STATS['total_reads_with_bc']:,}")
    print(f"Total reads corrected:            {STATS['total_corrected_reads']:,}")

    if STATS["total_reads_with_bc"] > 0:
        corr_pct = (STATS["total_corrected_reads"] / STATS["total_reads_with_bc"]) * 100
        print(f"Correction rate:                  {corr_pct:.2f}%")

    print("-" * 40)
    print(f"Unique barcodes (Before):         {len(STATS['unique_bcs_before']):,}")
    print(f"Unique barcodes (After):          {len(STATS['unique_bcs_after']):,}")

    unique_diff = len(STATS["unique_bcs_before"]) - len(STATS["unique_bcs_after"])
    print(f"Unique barcodes collapsed:        {unique_diff:,}")
    print("=" * 40 + "\n")


if __name__ == "__main__":
    main()
