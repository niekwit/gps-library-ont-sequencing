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
    counts = Counter(raw_bcs)
    unique_bcs = sorted(counts.keys(), key=lambda x: counts[x], reverse=True)

    # Update "Before" stats
    STATS["total_reads_with_bc"] += len(raw_bcs)
    STATS["unique_bcs_before"].update(unique_bcs)

    centroids = [b for b in unique_bcs if len(b) in args.target_lengths]
    correction_map = {}

    for seq in unique_bcs:
        if seq in centroids:
            correction_map[seq] = seq
            continue

        seq_count = counts[seq]
        candidates = []

        for c in centroids:
            if counts[c] >= (seq_count * args.ratio_threshold):
                if abs(len(seq) - len(c)) <= args.window:
                    dist = Levenshtein.distance(seq, c, score_cutoff=args.max_edit_dist)
                    if dist <= args.max_edit_dist:
                        candidates.append((dist, counts[c], c))

        if candidates:
            candidates.sort(key=lambda x: (x[0], -x[1]))
            correction_map[seq] = candidates[0][2]
        else:
            correction_map[seq] = seq

    # Write updated reads and update "After" stats
    for read in reads_with_bc:
        old_bc = read.get_tag("BC")
        new_bc = correction_map.get(old_bc, old_bc)

        if old_bc != new_bc:
            STATS["total_corrected_reads"] += 1
            read.set_tag("BC", new_bc)

        STATS["unique_bcs_after"].add(new_bc)
        writer.write(read)


def main():
    pysam.set_verbosity(0)

    with pysam.AlignmentFile(args.input, "rb") as reader:
        total_reads = reader.mapped
        with pysam.AlignmentFile(args.output, "wb", template=reader) as writer:
            current_rname = None
            reads_in_group = []

            for read in tqdm.tqdm(reader, total=total_reads, desc="Correcting Barcodes"):
                rname = read.reference_name
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
