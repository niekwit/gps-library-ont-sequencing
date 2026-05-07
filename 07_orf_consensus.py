"""
Strategy
--------
Reads in a barcode-sorted BAM of ONT reads aligned to Gencode cDNA sequences.
For each barcode (representing a library member), the script:

1. Groups reads by their BC tag and skips groups below MIN_READS.
2. Identifies the primary isoform — the reference transcript that the majority
   of reads in the group align to.
3. Computes a confidence score: the fraction of reads supporting that isoform.
4. Generates a pileup consensus from reads aligned to the primary isoform,
   recovering up to 10 bp of soft-clipped sequence on each end so that bases
   extending beyond the annotated transcript boundaries are not lost.
5. Fetches the matching Gencode cDNA sequence from the reference FASTA for
   direct comparison.
6. Writes one row per barcode to a CSV with the consensus, original cDNA,
   confidence score, and read counts.
"""

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


def extract_cds(header, seq):
    """Returns the CDS subsequence using coordinates from the Gencode header.

    Gencode headers contain a field like 'CDS:71-2188' with 1-based inclusive
    coordinates. Returns the full sequence unchanged if no CDS field is found.
    """
    for field in header.split("|"):
        if field.startswith("CDS:"):
            start, end = field[4:].split("-")
            # Convert from 1-based inclusive to 0-based Python slice
            return seq[int(start) - 1 : int(end)]
    return seq


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
        # Cache aligned pairs once per read to improve performance
        aligned_pairs = read.get_aligned_pairs()
        query_seq = read.query_sequence

        # Standard aligned bases
        for q_pos, r_pos in aligned_pairs:
            if r_pos is not None and 0 <= r_pos < ref_len:
                base = query_seq[q_pos] if q_pos is not None else "-"
                pileup[r_pos][base] += 1

        # Filter for pairs that actually map to the reference
        ref_mapped_indices = [
            i for i, p in enumerate(aligned_pairs) if p[1] is not None
        ]

        if not ref_mapped_indices:
            continue

        # Capture Prefixes
        first_idx = ref_mapped_indices[0]
        q_idx_start = aligned_pairs[first_idx][0]

        if q_idx_start is not None:
            for i in range(1, 11):
                if q_idx_start - i >= 0:
                    pre_base = query_seq[q_idx_start - i]
                    prefix_pileup[i - 1][pre_base] += 1

        # Capture Suffixes
        last_idx = ref_mapped_indices[-1]
        q_idx_end = aligned_pairs[last_idx][0]

        if q_idx_end is not None:
            for i in range(1, 11):
                if q_idx_end + i < len(query_seq):
                    suf_base = query_seq[q_idx_end + i]
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

            # Identify Majority Reference (Primary Isoform)
            ref_counts = Counter(r.reference_name for r in reads if r.reference_name)
            if not ref_counts:
                continue

            full_header, primary_count = ref_counts.most_common(1)[0]

            # Calculate Confidence Score
            confidence_score = round((primary_count / total_reads) * 100, 2)

            # Extract Gene/Transcript name (5th element of | separated header)
            header_parts = full_header.split("|")
            short_name = header_parts[5] if len(header_parts) >= 6 else "Unknown"

            # Fetch GenCode cDNA sequence and trim to CDS only
            try:
                original_cdna = extract_cds(full_header, ref_seqs.fetch(full_header))
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
