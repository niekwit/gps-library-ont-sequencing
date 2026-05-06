"""
Strategy
--------
For reads grouped by gene name (6th field of RNAME split by "|"):

1. Collect unique BC tag sequences from all reads in the group.
2. Compute pairwise Levenshtein distances across all barcode combinations.
3. If any pair falls within --threshold, those barcodes are flagged as
   "ambiguous" and enter the four-stage resolution pipeline below.

Resolution pipeline (each stage only runs if the previous produced no corrections):

  a. CSV lookup — fetch known barcodes for the gene from the reference CSV and
     run build_correction_map. Each ambiguous barcode is replaced by the nearest
     reference barcode within --threshold. Ties resolved by lowest distance;
     equal distances resolved by order of appearance in the CSV.

  b. Alias lookup — if the CSV lookup found no reference barcodes, query
     MyGene.info with the Ensembl gene ID (field 2 of RNAME, version stripped)
     to obtain alternative gene names. Each alias is checked against the
     reference CSV; the first hit is used as the reference set.

  c. Length-based fallback — if stages (a) and (b) produced no corrections,
     inspect the ambiguous barcodes themselves. If exactly one has a length of
     24 or 30 bp (the expected library barcode lengths), treat it as an internal
     anchor and attempt correction of the remaining ambiguous barcodes against it.

  d. Skip — if all three stages fail, log a WARNING and leave the reads unchanged.

4. All reads (corrected or not) are written to the output BAM. Corrected reads
   receive an updated BC tag, an OB tag preserving the original barcode, and
   XB=True. Uncorrected reads receive XB=False.
"""

import argparse
import logging
import mygene
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
# Logging — INFO goes to both console and file; all other levels to file only
# ---------------------------------------------------------------------------
_fmt = logging.Formatter(
    "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

_file_handler = logging.FileHandler(args.log_file)
_file_handler.setLevel(logging.DEBUG)
_file_handler.setFormatter(_fmt)


class _InfoOnly(logging.Filter):
    """Passes only INFO-level records — keeps WARNING and above off the console."""

    def filter(self, record):
        return record.levelno == logging.INFO


_console_handler = logging.StreamHandler()
_console_handler.setLevel(logging.INFO)
_console_handler.addFilter(_InfoOnly())
_console_handler.setFormatter(_fmt)

logging.basicConfig(level=logging.DEBUG, handlers=[_file_handler, _console_handler])
log = logging.getLogger(__name__)

# Silence the HTTP client libraries used internally by mygene so their
# connection/request details don't flood the log file.
logging.getLogger("httpcore").setLevel(logging.WARNING)
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("mygene").setLevel(logging.WARNING)

# ---------------------------------------------------------------------------
# Global metrics — accumulated across all gene groups and reported at the end
# ---------------------------------------------------------------------------
STATS = {
    "total_reads_with_bc": 0,
    "total_corrected_reads": 0,
    "unique_bcs_before": set(),  # all unique BC values seen before correction
    "unique_bcs_after": set(),  # all unique BC values written to the output
    "genes_skipped": 0,  # genes with ambiguous barcodes that could not be corrected
    "genes_alias_resolved": 0,  # genes resolved via Ensembl alias lookup
}

# ---------------------------------------------------------------------------
# Ensembl alias lookup — used when a gene name is absent from the reference CSV
# ---------------------------------------------------------------------------
_mg = mygene.MyGeneInfo()
_alias_cache: dict[str, list[str]] = {}  # ENSG ID (no version) -> aliases


def extract_ensg_id(rname):
    """Returns the versionless ENSG ID from a Gencode RNAME (2nd |-separated field)."""
    if not rname:
        return None
    parts = rname.split("|")
    if len(parts) < 2:
        return None
    # Strip the version suffix (e.g. ENSG00000179218.12 -> ENSG00000179218)
    return parts[1].split(".")[0]


def get_aliases(ensg_id):
    """Returns all known gene name aliases for an Ensembl gene ID via MyGene.info.

    Results are cached so each ENSG ID is only queried once per run.
    Returns an empty list if the query fails or the gene has no aliases.
    """
    if ensg_id in _alias_cache:
        return _alias_cache[ensg_id]

    try:
        result = _mg.getgene(
            ensg_id, fields=["alias", "symbol"], species="human", verbose=False
        )
        aliases = []
        if result:
            raw = result.get("alias", [])
            # mygene returns a str when there is only one alias, list otherwise
            aliases = [raw] if isinstance(raw, str) else list(raw)
            symbol = result.get("symbol")
            if symbol and symbol not in aliases:
                aliases.append(symbol)
    except Exception as e:
        log.debug(f"MyGene.info lookup failed for {ensg_id}: {e}")
        aliases = []

    _alias_cache[ensg_id] = aliases
    return aliases


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
            log.debug(f"Correction mapped: {bc} -> {best_bc} (distance {best_dist})")
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
            log.debug(f"Gene {gene}: {len(close_pairs)} close barcode pair(s) found")

            # Step 2: collect the full set of ambiguous barcodes (both members
            # of every close pair) — needed for both the reference lookup and
            # the length-based fallback below.
            ambiguous = set()
            for bc1, bc2, _ in close_pairs:
                ambiguous.add(bc1)
                ambiguous.add(bc2)

            # Step 3: look up the known ("ground truth") barcodes for this
            # gene from the reference CSV. Without a reference we have no
            # anchor to decide which of the two close sequences is correct.
            ref_bcs = ref_by_gene.get(gene, [])

            if ref_bcs:
                log.debug(f"Gene {gene}: found {len(ref_bcs)} reference barcode(s) in CSV")

            if not ref_bcs:
                # Alias fallback: the reference CSV may use an older name for
                # this gene. Query MyGene.info with the Ensembl gene ID and
                # check whether any alias matches an entry in the reference.
                ensg_id = extract_ensg_id(reads[0].reference_name)
                if ensg_id:
                    for alias in get_aliases(ensg_id):
                        ref_bcs = ref_by_gene.get(alias, [])
                        if ref_bcs:
                            log.debug(
                                f"Gene {gene}: found reference barcodes under "
                                f"alias '{alias}' ({ensg_id})"
                            )
                            STATS["genes_alias_resolved"] += 1
                            break

            if ref_bcs:
                correction_map = build_correction_map(
                    ambiguous, ref_bcs, args.threshold
                )
                if correction_map:
                    log.debug(f"Gene {gene}: {len(correction_map)} barcode(s) will be corrected")
                else:
                    log.debug(f"Gene {gene}: no ambiguous barcodes matched a reference within threshold")

            if not correction_map:
                # Length fallback: applies whenever CSV/alias refs either didn't
                # exist or existed but produced no corrections within threshold.
                # If exactly one ambiguous barcode has the expected length (24 or
                # 30), use it as an internal anchor to correct the others.
                target_len_bcs = [bc for bc in ambiguous if len(bc) in (24, 30)]
                if len(target_len_bcs) == 1:
                    log.debug(
                        f"Gene {gene}: using length-{len(target_len_bcs[0])} "
                        f"barcode {target_len_bcs[0]} as length-based fallback anchor."
                    )
                    correction_map = build_correction_map(
                        ambiguous, target_len_bcs, args.threshold
                    )
                    if correction_map:
                        log.debug(
                            f"Gene {gene}: {len(correction_map)} barcode(s) corrected "
                            f"via length-based fallback"
                        )
                    else:
                        log.debug(
                            f"Gene {gene}: length-based anchor produced no corrections within threshold"
                        )
                else:
                    STATS["genes_skipped"] += 1
                    log.warning(
                        f"Gene {gene} has {len(close_pairs)} ambiguous barcode "
                        f"pair(s) but no reference barcodes (CSV and alias lookup "
                        f"both failed) and no unambiguous length-24/30 anchor — "
                        f"skipping correction for this gene."
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
            read.set_tag("XB", True)  # flag: this read was corrected
        else:
            read.set_tag("XB", False)  # flag: no correction applied

        STATS["unique_bcs_after"].add(new_bc)
        writer.write(read)


def main():
    log.info("Starting reference barcode correction")
    log.info(f"Input BAM:          {args.input}")
    log.info(f"Output BAM:         {args.output}")
    log.info(f"Reference CSV:      {args.reference}")
    log.info(f"Distance threshold: {args.threshold}")
    log.info(f"Reverse complement: {args.rc}")
    log.info(f"Log file:           {args.log_file}")

    # ---------------------------------------------------------------------------
    # Load the reference barcode table and index it by gene name for O(1) lookup
    # inside process_group. Rows with missing gene or barcode values are dropped
    # because they cannot contribute to a meaningful correction.
    # ---------------------------------------------------------------------------
    log.info(f"Loading reference barcodes from: {args.reference}")
    ref_df = pd.read_csv(args.reference)
    ref_by_gene = (
        ref_df.dropna(subset=[args.gene_col, args.bc_col])
        .groupby(args.gene_col)[args.bc_col]
        .apply(list)
        .to_dict()
    )
    log.info(f"Loaded reference barcodes for {len(ref_by_gene)} genes.")

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
    log.info(f"Input file:            {args.input}")
    log.info(f"Output file:           {args.output}")
    log.info(f"Reference file:        {args.reference}")
    log.info(f"Distance threshold:    {args.threshold}")
    log.info(thin)
    log.info(f"Total reads (with BC): {STATS['total_reads_with_bc']:,}")
    log.info(f"Total reads corrected: {STATS['total_corrected_reads']:,}")

    if STATS["total_reads_with_bc"] > 0:
        pct = (STATS["total_corrected_reads"] / STATS["total_reads_with_bc"]) * 100
        log.info(f"Correction rate:       {pct:.2f}%")

    log.info(thin)
    log.info(f"Unique BCs (before):   {len(STATS['unique_bcs_before']):,}")
    log.info(f"Unique BCs (after):    {len(STATS['unique_bcs_after']):,}")
    collapsed = len(STATS["unique_bcs_before"]) - len(STATS["unique_bcs_after"])
    log.info(f"Unique BCs collapsed:  {collapsed:,}")
    log.info(thin)
    log.info(
        f"Genes resolved via alias lookup:          {STATS['genes_alias_resolved']:,}"
    )
    log.info(f"Genes skipped (ambiguous, uncorrectable): {STATS['genes_skipped']:,}")
    log.info(divider)


if __name__ == "__main__":
    main()
