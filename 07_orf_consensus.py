import argparse
import pandas as pd
import tqdm
import pysam
from itertools import groupby
from collections import Counter

# Command line parser
# ----------------------------
parser = argparse.ArgumentParser(
    description="Generate ORF consensus sequences with confidence scoring"
)
parser.add_argument(
    "-i", "--input", required=True, help="Input BAM file (sorted by BC tag)"
)
parser.add_argument(
    "-o", "--output", help="Output CSV file", default="raw_consensus_sequences.csv"
)
parser.add_argument(
    "-m",
    "--min-reads",
    type=int,
    help="Minimum reads per barcode to generate a consensus",
    default=3,
)
parser.add_argument(
    "-f", "--ref-fasta", required=True, help="Gencode Reference FASTA file"
)

# Configuration
# ----------------------------
args = parser.parse_args()
BAM_FILE = args.input
REF_FASTA = args.ref_fasta
MIN_READS = args.min_reads
OUTPUT_CSV = args.output


def get_bc_tag(read):
    """Retrieves the BC tag from a pysam AlignedSegment."""
    return read.get_tag("BC") if read.has_tag("BC") else None


def get_pileup_consensus(reads, ref_len):
    """Generates consensus sequence including soft-clipped starts and ends."""
    pileup = [Counter() for _ in range(ref_len)]

    # Counters for bases appearing before/after the reference range
    prefix_pileup = [Counter() for _ in range(10)]  # Check up to 10bp overhang
    suffix_pileup = [Counter() for _ in range(10)]

    for read in reads:
        # Standard aligned bases
        for q_pos, r_pos in read.get_aligned_pairs():
            if r_pos is not None and 0 <= r_pos < ref_len:
                base = read.query_sequence[q_pos] if q_pos is not None else "-"
                pileup[r_pos][base] += 1

        # Capture Prefixes (a potentially "lost" ATG)
        # Find the first aligned reference position for this read
        aligned_positions = [p[1] for p in read.get_aligned_pairs() if p[1] is not None]
        if aligned_positions:
            first_ref_idx = read.get_aligned_pairs().index(
                (
                    next(
                        p[0]
                        for p in read.get_aligned_pairs()
                        if p[1] == aligned_positions[0]
                    ),
                    aligned_positions[0],
                )
            )

            # Look backwards from the first aligned base in the read
            for i in range(1, 11):
                q_idx = read.get_aligned_pairs()[first_ref_idx][0]
                if q_idx is not None and q_idx - i >= 0:
                    pre_base = read.query_sequence[q_idx - i]
                    prefix_pileup[i - 1][pre_base] += 1

        # Capture Suffixes (the "lost" Stop)
        if aligned_positions:
            last_ref_idx = read.get_aligned_pairs().index(
                (
                    next(
                        p[0]
                        for p in reversed(read.get_aligned_pairs())
                        if p[1] == aligned_positions[-1]
                    ),
                    aligned_positions[-1],
                )
            )

            for i in range(1, 11):
                q_idx = read.get_aligned_pairs()[last_ref_idx][0]
                if q_idx is not None and q_idx + i < len(read.query_sequence):
                    suf_base = read.query_sequence[q_idx + i]
                    suffix_pileup[i - 1][suf_base] += 1

    # Build the final sequence
    res = []

    # Add prefix (reversed because we counted backwards from the start)
    for counts in reversed(prefix_pileup):
        if counts and counts.most_common(1)[0][1] >= (len(reads) / 2):
            res.append(counts.most_common(1)[0][0])

    # Add middle (standard pileup)
    for counts in pileup:
        if counts:
            best = counts.most_common(1)[0][0]
            if best != "-":
                res.append(best)

    # Add suffix
    for counts in suffix_pileup:
        if counts and counts.most_common(1)[0][1] >= (len(reads) / 2):
            res.append(counts.most_common(1)[0][0])

    return "".join(res)


def process_to_csv():
    # Load reference FASTA for the "original" sequences
    print(f"Loading reference FASTA: {REF_FASTA}")
    ref_seqs = pysam.FastaFile(REF_FASTA)

    results = []

    with pysam.AlignmentFile(BAM_FILE, "rb", check_sq=False) as bam:
        # fetch(until_eof=True) allows streaming tag-sorted BAMs without an index
        reader = bam.fetch(until_eof=True)

        print("Processing barcodes and generating consensus...")
        # Use tqdm to track progress through the iterator
        for bc, group in tqdm.tqdm(groupby(reader, key=get_bc_tag)):
            if not bc:
                continue

            # Materialize group into a list to allow multiple passes
            reads = list(group)
            total_reads = len(reads)

            if total_reads < MIN_READS:
                continue

            # 1. Identify Majority Reference (Primary Isoform)
            ref_counts = Counter(r.reference_name for r in reads if r.reference_name)
            if not ref_counts:
                continue

            full_header, primary_count = ref_counts.most_common(1)[0]

            # Calculate Confidence Score
            confidence_score = round((primary_count / total_reads) * 100, 2)

            # Extract Gene/Transcript name (5th element of | separated header)
            header_parts = full_header.split("|")
            short_name = header_parts[4] if len(header_parts) >= 5 else "Unknown"

            # Fetch GenCode cDNA sequence (from FASTA)
            try:
                original_cdna = ref_seqs.fetch(full_header)
            except KeyError:
                original_cdna = "REF_NOT_FOUND"

            # Generate consensus using ONLY reads aligned to the primary reference
            ref_len = bam.header.get_reference_length(full_header)
            primary_reads = [r for r in reads if r.reference_name == full_header]
            consensus = get_pileup_consensus(primary_reads, ref_len)

            results.append(
                {
                    "barcode": bc,
                    "gene_name": short_name,
                    "gencode_header": full_header,
                    "total_reads": total_reads,
                    "primary_reads": primary_count,
                    "confidence_pct": confidence_score,
                    "original_cdna": original_cdna,
                    "consensus_cdna": consensus,
                }
            )

    # Save to CSV
    print(f"Aggregating {len(results)} results into DataFrame...")
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Done! Results saved to {OUTPUT_CSV}")


if __name__ == "__main__":
    process_to_csv()
