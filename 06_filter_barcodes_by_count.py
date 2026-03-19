import pysam
import sys
from collections import Counter
from tqdm import tqdm


def filter_bam_by_bc_count(input_bam, output_bam, min_count=5):
    # Pass 1: Count unique barcodes
    bc_counts = Counter()

    print(f"Pass 1: Counting barcodes in {input_bam}...")
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        # Get total reads for tqdm progress bar
        total_reads = bam.mapped + bam.unmapped
        for read in tqdm(bam, total=total_reads, desc="Counting"):
            if read.has_tag("BC"):
                bc = read.get_tag("BC")
                bc_counts[bc] += 1

    # Pass 2: Filter and write to new BAM
    print(f"Pass 2: Filtering reads (min_count={min_count})...")
    with pysam.AlignmentFile(input_bam, "rb") as reader:
        # Use 'template=reader' to copy the header (references, etc.)
        with pysam.AlignmentFile(output_bam, "wb", template=reader) as writer:
            for read in tqdm(reader, total=total_reads, desc="Filtering"):
                if read.has_tag("BC"):
                    bc = read.get_tag("BC")
                    if bc_counts[bc] >= min_count:
                        writer.write(read)

    print(f"Done! Filtered BAM saved to: {output_bam}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 07_filter_bam.py <input.bam> <output.bam> [min_count]")
        sys.exit(1)

    in_bam = sys.argv[1]
    out_bam = sys.argv[2]
    threshold = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    filter_bam_by_bc_count(in_bam, out_bam, threshold)
