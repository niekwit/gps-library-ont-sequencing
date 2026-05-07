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
   e. If multiple ground-truth candidates are found, use read support as a
      tiebreaker: the candidate with the most rows in the gene group is
      chosen.  If two candidates share the exact same row count the
      component is skipped and a warning is logged.
   f. If zero ground-truth candidates are found, skip the component and log
      a warning.
5. Write the corrected DataFrame to a new CSV with an added orf_corrected
   boolean column.
"""

import argparse
import logging
from collections import defaultdict

import pandas as pd
import tqdm
from rapidfuzz.distance import Levenshtein

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Correct ORF consensus sequences using Levenshtein distance and ORF validity."
)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input CSV file (output of 07_orf_consensus.py)",
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
    default=2,
    help="Maximum Levenshtein distance to consider two sequences as related (default: 2)",
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
    """Returns all (seq_a, seq_b, distance) pairs with distance ≤ threshold.

    Sequences are sorted by length so the inner loop can break early: once the
    length difference alone exceeds the threshold, Levenshtein distance is
    guaranteed to exceed it too, so all remaining comparisons can be skipped.
    """
    pairs = []
    sorted_seqs = sorted(seqs, key=len)
    n = len(sorted_seqs)
    for i in range(n):
        for j in range(i + 1, n):
            if len(sorted_seqs[j]) - len(sorted_seqs[i]) > threshold:
                break
            dist = Levenshtein.distance(
                sorted_seqs[i], sorted_seqs[j], score_cutoff=threshold
            )
            if dist <= threshold:
                pairs.append((sorted_seqs[i], sorted_seqs[j], dist))
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
        "components_multi_truth_resolved": 0,
        "rows_corrected": 0,
    }

    # Collect all corrections during the read-only loop pass; apply after.
    # Modifying df inside a groupby iteration triggers pandas CoW invalidation
    # on every write, which causes O(n_genes * n_rows) overhead.
    all_corrections: dict[int, str] = {}  # DataFrame row index → corrected sequence

    n_genes = df["gene_name"].nunique()
    for gene_name, group in tqdm.tqdm(
        df.groupby("gene_name"), total=n_genes, desc="Processing genes"
    ):
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

        correction_map: dict[str, str] = {}

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
                # Tiebreaker: pick the valid ORF with the most rows in this group
                seq_counts = group["consensus_cdna"].value_counts()
                valid_orfs_sorted = sorted(
                    valid_orfs, key=lambda s: seq_counts.get(s, 0), reverse=True
                )
                top = seq_counts.get(valid_orfs_sorted[0], 0)
                second = seq_counts.get(valid_orfs_sorted[1], 0)
                if top == second:
                    stats["components_multi_truth"] += 1
                    log.warning(
                        f"Gene {gene_name}: {len(valid_orfs)} valid ORF candidates tied "
                        f"on barcode count — skipping. Candidates: {valid_orfs}"
                    )
                    continue
                stats["components_multi_truth_resolved"] += 1
                gt_preview = valid_orfs_sorted[0]
                log.debug(
                    f"Gene {gene_name}: {len(valid_orfs)} valid ORF candidates resolved "
                    f"by barcode count ({top} vs {second}) — "
                    f"chose '{gt_preview[:40]}{'...' if len(gt_preview) > 40 else ''}'"
                )
                valid_orfs = [valid_orfs_sorted[0]]

            ground_truth = valid_orfs[0]
            for seq in comp - {ground_truth}:
                correction_map[seq] = ground_truth

            log.debug(
                f"Gene {gene_name}: correcting {len(comp) - 1} sequence(s) → "
                f"ground truth '{ground_truth[:40]}{'...' if len(ground_truth) > 40 else ''}'"
            )
            stats["components_corrected"] += 1

        if correction_map:
            for row_idx, old_seq in group["consensus_cdna"].items():
                if old_seq in correction_map:
                    all_corrections[row_idx] = correction_map[old_seq]

    # Apply all corrections in one batch — no per-gene DataFrame writes
    if all_corrections:
        idx = list(all_corrections.keys())
        new_seqs = [all_corrections[i] for i in idx]
        df.loc[idx, "consensus_cdna"] = new_seqs
        df.loc[idx, "consensus_length"] = [len(s) for s in new_seqs]
        df.loc[idx, "orf_corrected"] = True
        stats["rows_corrected"] = len(idx)

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
    log.info(f"ORF clusters found:                         {stats['components_found']:,}")
    log.info(f"ORF clusters corrected:                     {stats['components_corrected']:,}")
    log.info(f"ORF clusters — no valid ORF:                {stats['components_no_truth']:,}")
    log.info(f"ORF clusters — ambiguous (resolved by BC):  {stats['components_multi_truth_resolved']:,}")
    log.info(f"ORF clusters — ambiguous (tied, skipped):   {stats['components_multi_truth']:,}")
    log.info(divider)


if __name__ == "__main__":
    process()
