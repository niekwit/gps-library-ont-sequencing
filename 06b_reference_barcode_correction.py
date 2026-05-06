"""
Strategy
--------
For reads grouped by gene name (6th field of RNAME split by "|"):

1. Collect unique BC tag sequences from all reads in the group.
2. Compute pairwise Levenshtein distances across all barcode combinations.
3. If any pair falls within --threshold, those barcodes are "ambiguous" and
   need adjudication against a reference barcode set.
4. Fetch the reference barcodes for that gene from the supplied CSV.
5. For each ambiguous barcode, find the nearest reference barcode. If one is
   within --threshold, update the BC tag on every affected read to that
   reference sequence. Ties are broken by choosing the lower distance; if
   distances are equal the first reference barcode wins.
6. All reads (corrected or not) are written to the output BAM. Corrected reads
   also receive an OB tag (original barcode) and XB flag set to True.
"""

import argparse
import logging
import pysam
import pandas as pd
from itertools import combinations
from rapidfuzz.distance import Levenshtein
from Bio.Seq import Seq
import tqdm

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Correct ambiguous barcodes against a reference barcode set."
)
parser.add_argument(
    "-i", "--input", required=True, help="Input BAM file (sorted by RNAME)"
)
parser.add_argument("-o", "--output", required=True, help="Output BAM file")
parser.add_argument(
    "-r",
    "--reference",
    required=True,
    help="Reference CSV file containing gene names and known barcode sequences",
)
parser.add_argument(
    "--gene-col",
    help="Column name for gene names in the reference CSV (default: gene_name_x)",
)
parser.add_argument(
    "--bc-col",
    help="Column name for barcode sequences in the reference CSV (default: barcode)",
)
parser.add_argument(
    "--threshold",
    type=int,
    default=2,
    help="Levenshtein distance threshold to flag a barcode pair as ambiguous (default: 2)",
)
parser.add_argument(
    "--rc",
    action="store_true",
    help="Reverse complement all reference barcode sequences after loading",
)
parser.add_argument(
    "--log-file",
    default="barcode_correction.log",
    help="Path to the log file (default: barcode_correction.log)",
)

args = parser.parse_args()

# ---------------------------------------------------------------------------
# Logging — INFO and above go to both console and file; DEBUG goes to file only
# ---------------------------------------------------------------------------
_fmt = logging.Formatter(
    "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

_file_handler = logging.FileHandler(args.log_file)
_file_handler.setLevel(logging.DEBUG)
_file_handler.setFormatter(_fmt)

_console_handler = logging.StreamHandler()
_console_handler.setLevel(logging.INFO)
_console_handler.setFormatter(_fmt)

logging.basicConfig(level=logging.DEBUG, handlers=[_file_handler, _console_handler])
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Global metrics — accumulated across all gene groups and reported at the end
# ---------------------------------------------------------------------------
STATS = {
    "total_reads_with_bc": 0,
    "total_corrected_reads": 0,
    "unique_bcs_before": set(),  # all unique BC values seen before correction
    "unique_bcs_after": set(),   # all unique BC values written to the output
}


def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence using Biopython."""
    return str(Seq(seq).reverse_complement())


def extract_gene_name(rname):
    """Returns the gene name from a Gencode RNAME (6th |-separated field)."""
    if not rname:
        return None
    parts = rname.split("|")
    # Gencode headers follow the pattern:
    # ENST...|ENSG...|OTTHUM...|OTTHUM...|transcript_name|gene_name|...
    return parts[5] if len(parts) >= 6 else rname


def find_close_pairs(barcodes, threshold):
    """Returns (bc1, bc2, distance) tuples for all pairs within threshold.

    Uses score_cutoff so rapidfuzz can exit early once a partial alignment
    exceeds the threshold, avoiding unnecessary work on clearly distant pairs.
    """
    result = []
    for bc1, bc2 in combinations(barcodes, 2):
        d = Levenshtein.distance(bc1, bc2, score_cutoff=threshold)
        if d <= threshold:
            result.append((bc1, bc2, d))
    return result


def build_correction_map(ambiguous_bcs, ref_bcs, threshold):
    """Maps each ambiguous barcode to its nearest reference barcode.

    Only barcodes that actually change (old != new) are included in the
    returned dict, so the caller can use .get(bc, bc) as a no-op default.

    Ties (equal distance to two reference barcodes) resolve to whichever
    reference barcode appears first in ref_bcs — deterministic but arbitrary.
    """
    correction_map = {}
    for bc in ambiguous_bcs:
        best_bc = None
        best_dist = threshold + 1  # sentinel: one above the allowed maximum

        for ref_bc in ref_bcs:
            d = Levenshtein.distance(bc, ref_bc, score_cutoff=threshold)
            # Keep only if strictly better than the current best
            if d <= threshold and d < best_dist:
                best_dist = d
                best_bc = ref_bc

        # Only record a correction when a better sequence was actually found
        if best_bc is not None and best_bc != bc:
            log.debug("Correction mapped: %s -> %s (distance %d)", bc, best_bc, best_dist)
            correction_map[bc] = best_bc

    return correction_map


def process_group(reads, writer, ref_by_gene):
    """Corrects ambiguous barcodes for one gene group and writes all reads."""

    # Partition reads: those without a BC tag cannot be corrected and are
    # passed straight through to the output unchanged.
    reads_with_bc = [r for r in reads if r.has_tag("BC")]
    for r in reads:
        if not r.has_tag("BC"):
            writer.write(r)

    if not reads_with_bc:
        return

    # Collect raw BC strings and derive the set of unique sequences seen in
    # this gene group. Deduplication is needed because we compare barcodes
    # against each other, not individual reads.
    raw_bcs = [r.get_tag("BC") for r in reads_with_bc]
    unique_bcs = list(set(raw_bcs))

    # Update global "before" metrics before any corrections are applied
    STATS["total_reads_with_bc"] += len(raw_bcs)
    STATS["unique_bcs_before"].update(unique_bcs)

    correction_map = {}  # populated below only if ambiguous pairs exist

    if len(unique_bcs) >= 2:
        # Step 1: find all barcode pairs that are suspiciously similar.
        # A pair within the threshold could be a sequencing error — one of
        # them might be a corrupted version of the true library barcode.
        close_pairs = find_close_pairs(unique_bcs, args.threshold)

        if close_pairs:
            gene = extract_gene_name(reads[0].reference_name)
            log.debug("Gene %s: %d close barcode pair(s) found", gene, len(close_pairs))

            # Step 2: look up the known ("ground truth") barcodes for this
            # gene from the reference CSV. Without a reference we have no
            # anchor to decide which of the two close sequences is correct.
            ref_bcs = ref_by_gene.get(gene, [])

            if not ref_bcs:
                log.warning(
                    "Gene %s has %d ambiguous barcode pair(s) but no reference "
                    "barcodes — skipping correction for this gene.",
                    gene,
                    len(close_pairs),
                )
            else:
                # Step 3: collect the full set of ambiguous barcodes (both
                # members of every close pair) and find the nearest reference
                # barcode for each one.
                ambiguous = set()
                for bc1, bc2, _ in close_pairs:
                    ambiguous.add(bc1)
                    ambiguous.add(bc2)

                correction_map = build_correction_map(
                    ambiguous, ref_bcs, args.threshold
                )

    # Step 4: apply corrections and write every read to the output BAM.
    for read in reads_with_bc:
        old_bc = read.get_tag("BC")
        # Falls back to old_bc unchanged if no correction was found
        new_bc = correction_map.get(old_bc, old_bc)

        # Always preserve the original barcode for traceability
        read.set_tag("OB", old_bc)

        if old_bc != new_bc:
            STATS["total_corrected_reads"] += 1
            read.set_tag("BC", new_bc)
            read.set_tag("XB", True)   # flag: this read was corrected
        else:
            read.set_tag("XB", False)  # flag: no correction applied

        STATS["unique_bcs_after"].add(new_bc)
        writer.write(read)


def main():
    log.info("Starting reference barcode correction")
    log.info("Input BAM:          %s", args.input)
    log.info("Output BAM:         %s", args.output)
    log.info("Reference CSV:      %s", args.reference)
    log.info("Distance threshold: %d", args.threshold)
    log.info("Reverse complement: %s", args.rc)
    log.info("Log file:           %s", args.log_file)

    # ---------------------------------------------------------------------------
    # Load the reference barcode table and index it by gene name for O(1) lookup
    # inside process_group. Rows with missing gene or barcode values are dropped
    # because they cannot contribute to a meaningful correction.
    # ---------------------------------------------------------------------------
    log.info("Loading reference barcodes from: %s", args.reference)
    ref_df = pd.read_csv(args.reference)
    ref_by_gene = (
        ref_df.dropna(subset=[args.gene_col, args.bc_col])
        .groupby(args.gene_col)[args.bc_col]
        .apply(list)
        .to_dict()
    )
    log.info("Loaded reference barcodes for %d genes.", len(ref_by_gene))

    # Reverse complement all reference barcodes when the BC sequences in the
    # BAM were captured from the opposite strand to the reference CSV.
    if args.rc:
        ref_by_gene = {
            gene: [reverse_complement(bc) for bc in bcs]
            for gene, bcs in ref_by_gene.items()
        }
        log.info("Reverse complemented all reference barcode sequences.")

    pysam.set_verbosity(0)

    with pysam.AlignmentFile(args.input, "rb") as reader:
        with pysam.AlignmentFile(args.output, "wb", template=reader) as writer:

            # ---------------------------------------------------------------------------
            # Stream through the BAM one read at a time and accumulate reads into a
            # group as long as they share the same gene name. When the gene changes,
            # flush the completed group for processing before starting the new one.
            # This avoids loading the entire BAM into memory while still enabling
            # within-gene barcode comparisons.
            # The input BAM must be sorted by RNAME (i.e. all reads for one gene are
            # consecutive) — this is a natural consequence of aligning to cDNA refs.
            # ---------------------------------------------------------------------------
            current_gene = None
            reads_in_group = []

            for read in tqdm.tqdm(reader, desc="Processing reads"):
                gene = extract_gene_name(read.reference_name)

                if gene != current_gene:
                    # Gene boundary reached — process the completed group
                    if reads_in_group:
                        process_group(reads_in_group, writer, ref_by_gene)
                    current_gene = gene
                    reads_in_group = [read]
                else:
                    reads_in_group.append(read)

            # Flush the final group which never triggers the boundary check above
            if reads_in_group:
                process_group(reads_in_group, writer, ref_by_gene)

    # ---------------------------------------------------------------------------
    # Summary report
    # ---------------------------------------------------------------------------
    divider = "=" * 40
    thin = "-" * 40
    log.info(divider)
    log.info("REFERENCE BARCODE CORRECTION SUMMARY")
    log.info(divider)
    log.info("Input file:            %s", args.input)
    log.info("Output file:           %s", args.output)
    log.info("Reference file:        %s", args.reference)
    log.info("Distance threshold:    %d", args.threshold)
    log.info(thin)
    log.info("Total reads (with BC): %d", STATS["total_reads_with_bc"])
    log.info("Total reads corrected: %d", STATS["total_corrected_reads"])

    if STATS["total_reads_with_bc"] > 0:
        pct = (STATS["total_corrected_reads"] / STATS["total_reads_with_bc"]) * 100
        log.info("Correction rate:       %.2f%%", pct)

    log.info(thin)
    log.info("Unique BCs (before):   %d", len(STATS["unique_bcs_before"]))
    log.info("Unique BCs (after):    %d", len(STATS["unique_bcs_after"]))
    collapsed = len(STATS["unique_bcs_before"]) - len(STATS["unique_bcs_after"])
    log.info("Unique BCs collapsed:  %d", collapsed)
    log.info(divider)


if __name__ == "__main__":
    main()
