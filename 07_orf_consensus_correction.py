"""
Post-process the consensus ORF CSV produced by 07_orf_consensus.py.

Strategy
--------
1. Load the consensus CSV (output of 07_orf_consensus.py, optionally extended
   with orf_id / barcode_id / orf_library_id / cds_ref columns).
2. Append the ORF8.1 linker sequence (CCAACTTTCTTGTACAAAGTGGTTTAA) to the
   consensus_cdna of every row whose orf_library_id is "ORF8.1".
3. Add a consensus_length column (length of consensus_cdna after any linker
   has been appended).
4. Group rows by gene_name.  Within each gene group:
   a. Find all pairs of distinct consensus_cdna sequences whose Levenshtein
      distance is ≤ --threshold.
   b. Build a graph where each node is a unique consensus sequence and each
      edge connects two sequences within the threshold.  Identify connected
      components via BFS.
   c. For each component with ≥ 2 members, attempt to find a single "ground
      truth" sequence defined as:
        - length that is a multiple of 3
        - starts with ATG
        - ends with a stop codon (TAA, TAG, or TGA)
   d. If exactly one ground-truth candidate exists within the component,
      replace every other sequence in that component with the ground truth
      and mark those rows as corrected (orf_corrected = True).
   e. If zero or multiple ground-truth candidates are found, skip the
      component and log a warning.
5. Write the corrected DataFrame to a new CSV with an added orf_corrected
   boolean column.
"""

import argparse
import logging
from collections import defaultdict
from itertools import combinations

import pandas as pd
from rapidfuzz.distance import Levenshtein

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Correct ORF consensus sequences using Levenshtein distance and ORF validity."
)
parser.add_argument(
    "-i", "--input", required=True, help="Input CSV file (output of 07_orf_consensus.py)"
)
parser.add_argument(
    "-o",
    "--output",
    default="orf_consensus_corrected.csv",
    help="Output CSV file (default: orf_consensus_corrected.csv)",
)
parser.add_argument(
    "-t",
    "--threshold",
    type=int,
    default=5,
    help="Maximum Levenshtein distance to consider two sequences as related (default: 5)",
)
parser.add_argument(
    "--log-file",
    default="orf_consensus_correction.log",
    help="Path to the log file (default: orf_consensus_correction.log)",
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

LINKER_ORF81 = "CCAACTTTCTTGTACAAAGTGGTTTAA"
STOP_CODONS = {"TAA", "TAG", "TGA"}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def is_valid_orf(seq: str) -> bool:
    """Returns True if seq is a plausible ground-truth ORF.

    Criteria: length divisible by 3, starts with ATG, ends with a stop codon.
    """
    if len(seq) % 3 != 0:
        return False
    if not seq.startswith("ATG"):
        return False
    if seq[-3:] not in STOP_CODONS:
        return False
    return True


def find_close_pairs(seqs: list[str], threshold: int) -> list[tuple[str, str, int]]:
    """Returns all (seq_a, seq_b, distance) pairs with distance ≤ threshold."""
    pairs = []
    for a, b in combinations(seqs, 2):
        dist = Levenshtein.distance(a, b, score_cutoff=threshold)
        if dist <= threshold:
            pairs.append((a, b, dist))
    return pairs


def connected_components(pairs: list[tuple[str, str, int]]) -> list[set[str]]:
    """Returns connected components from a list of (node1, node2, weight) edges."""
    graph: dict[str, set[str]] = defaultdict(set)
    for a, b, _ in pairs:
        graph[a].add(b)
        graph[b].add(a)

    visited: set[str] = set()
    components: list[set[str]] = []

    for start in graph:
        if start in visited:
            continue
        component: set[str] = set()
        queue = [start]
        while queue:
            node = queue.pop()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            queue.extend(graph[node] - visited)
        components.append(component)

    return components


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------


def process():
    log.info(f"Input CSV:   {args.input}")
    log.info(f"Output CSV:  {args.output}")
    log.info(f"Threshold:   {args.threshold}")
    log.info(f"Log file:    {args.log_file}")

    df = pd.read_csv(args.input)
    log.info(f"Loaded {len(df):,} rows from {args.input}")

    # Append linker for ORF8.1 entries before calculating lengths
    if "orf_library_id" in df.columns:
        mask_81 = df["orf_library_id"] == "ORF8.1"
        n_81 = mask_81.sum()
        if n_81 > 0:
            df.loc[mask_81, "consensus_cdna"] = (
                df.loc[mask_81, "consensus_cdna"] + LINKER_ORF81
            )
            log.info(f"Appended ORF8.1 linker to {n_81:,} rows")

    df["consensus_length"] = df["consensus_cdna"].str.len()

    # Initialise correction tracking columns
    df["orf_corrected"] = False

    # ------------------------------------------------------------------
    # Process each gene group
    # ------------------------------------------------------------------
    stats = {
        "components_found": 0,
        "components_corrected": 0,
        "components_no_truth": 0,
        "components_multi_truth": 0,
        "rows_corrected": 0,
    }

    for gene_name, group in df.groupby("gene_name"):
        unique_seqs = group["consensus_cdna"].dropna().unique().tolist()

        if len(unique_seqs) < 2:
            # Nothing to compare within a single-sequence gene
            continue

        pairs = find_close_pairs(unique_seqs, args.threshold)
        if not pairs:
            log.debug(f"Gene {gene_name}: no sequences within threshold, skipping")
            continue

        components = connected_components(pairs)
        log.debug(
            f"Gene {gene_name}: {len(unique_seqs)} unique sequences, "
            f"{len(pairs)} close pair(s), {len(components)} component(s)"
        )
        stats["components_found"] += len(components)

        for comp in components:
            if len(comp) < 2:
                continue

            valid_orfs = [s for s in comp if is_valid_orf(s)]

            if len(valid_orfs) == 0:
                stats["components_no_truth"] += 1
                log.warning(
                    f"Gene {gene_name}: component with {len(comp)} sequence(s) has no "
                    f"valid ORF — skipping. Sequences: {sorted(comp)}"
                )
                continue

            if len(valid_orfs) > 1:
                stats["components_multi_truth"] += 1
                log.warning(
                    f"Gene {gene_name}: component with {len(comp)} sequence(s) has "
                    f"{len(valid_orfs)} valid ORF candidates — ambiguous, skipping. "
                    f"Candidates: {valid_orfs}"
                )
                continue

            ground_truth = valid_orfs[0]
            to_correct = comp - {ground_truth}

            log.debug(
                f"Gene {gene_name}: correcting {len(to_correct)} sequence(s) → "
                f"ground truth '{ground_truth[:40]}{'...' if len(ground_truth) > 40 else ''}'"
            )

            # Apply correction to all matching rows in this gene's group
            mask = (df["gene_name"] == gene_name) & (df["consensus_cdna"].isin(to_correct))
            n_corrected = mask.sum()
            df.loc[mask, "consensus_cdna"] = ground_truth
            df.loc[mask, "consensus_length"] = len(ground_truth)
            df.loc[mask, "orf_corrected"] = True

            stats["components_corrected"] += 1
            stats["rows_corrected"] += n_corrected

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    df.to_csv(args.output, index=False)

    divider = "=" * 40
    thin = "-" * 40
    log.info(divider)
    log.info("ORF CONSENSUS CORRECTION SUMMARY")
    log.info(divider)
    log.info(f"Input CSV:                    {args.input}")
    log.info(f"Output CSV:                   {args.output}")
    log.info(f"Levenshtein threshold:        {args.threshold}")
    log.info(thin)
    log.info(f"Total rows:                   {len(df):,}")
    log.info(f"Rows corrected:               {stats['rows_corrected']:,}")
    log.info(thin)
    log.info(f"Components found:             {stats['components_found']:,}")
    log.info(f"Components corrected:         {stats['components_corrected']:,}")
    log.info(f"Components — no valid ORF:    {stats['components_no_truth']:,}")
    log.info(f"Components — ambiguous ORF:   {stats['components_multi_truth']:,}")
    log.info(divider)


if __name__ == "__main__":
    process()
