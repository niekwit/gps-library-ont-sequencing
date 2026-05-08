"""
Match consensus ORF sequences against GPS reference library entries.

Strategy
--------
1. Load both reference libraries:
     - ORFeome81 (annotations/ORFeome81_library_withAnnotation.csv):
       gene keyed by Gene_Symbol, IDs in ORF_ID
     - uORF (annotations/uORF_library_withAnnotation.csv):
       gene keyed by Gene Symbol, IDs in Clone ID
2. Read the consensus CSV produced by 07_orf_consensus_correction.py.
3. Group rows by gene_name.  Within each gene group:
   a. Determine which library to query from the most common barcode length:
        30 ± 2 bp → ORFeome81
        24 ± 2 bp → uORF
   b. Look up the gene name in the chosen library.  If not found, query
      MyGene.info for aliases and retry.
   c. For each unique original_cdna value in the group, find the reference
      entry whose nt_seq is within --threshold Levenshtein distance.  The
      closest match is chosen; ties are broken by higher ORF_ID/Clone_ID.
4. Write the matched reference ID to a new column matched_ref_id (NA when
   no match is found within the threshold).
"""

import argparse
import logging

import mygene
import pandas as pd
import tqdm
from rapidfuzz.distance import Levenshtein

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Match consensus ORF sequences against GPS reference libraries."
)
parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input CSV (output of 07_orf_consensus_correction.py)",
)
parser.add_argument(
    "-o",
    "--output",
    default="orfs_with_ref_ids.csv",
    help="Output CSV file (default: orfs_with_ref_ids.csv)",
)
parser.add_argument(
    "--orf81-ref",
    default="annotations/ORFeome81_library_withAnnotation.csv",
    help="ORFeome81 reference CSV (default: annotations/ORFeome81_library_withAnnotation.csv)",
)
parser.add_argument(
    "--uorf-ref",
    default="annotations/uORF_library_withAnnotation.csv",
    help="uORF reference CSV (default: annotations/uORF_library_withAnnotation.csv)",
)
parser.add_argument(
    "-t",
    "--threshold",
    type=int,
    default=5,
    help="Maximum Levenshtein distance for a sequence match (default: 5)",
)
parser.add_argument(
    "--log-file",
    default="reference_orf_matching.log",
    help="Path to the log file (default: reference_orf_matching.log)",
)

args = parser.parse_args()

# ---------------------------------------------------------------------------
# Logging — INFO to console + file; DEBUG/WARNING to file only
# ---------------------------------------------------------------------------
_fmt = logging.Formatter(
    "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

_file_handler = logging.FileHandler(args.log_file)
_file_handler.setLevel(logging.DEBUG)
_file_handler.setFormatter(_fmt)


class _InfoOnly(logging.Filter):
    def filter(self, record):
        return record.levelno == logging.INFO


_console_handler = logging.StreamHandler()
_console_handler.setLevel(logging.INFO)
_console_handler.addFilter(_InfoOnly())
_console_handler.setFormatter(_fmt)

logging.basicConfig(level=logging.DEBUG, handlers=[_file_handler, _console_handler])
log = logging.getLogger(__name__)

# Suppress noisy HTTP loggers from mygene / httpx
for _noisy in ("httpcore", "httpx", "mygene"):
    logging.getLogger(_noisy).setLevel(logging.WARNING)

# ---------------------------------------------------------------------------
# Load reference libraries
# ---------------------------------------------------------------------------


def load_references(orf81_path: str, uorf_path: str):
    """Returns two dicts: gene → list of (nt_seq, ref_id).

    orf81_ref : keyed by Gene_Symbol, ref_id = ORF_ID
    uorf_ref  : keyed by Gene Symbol, ref_id = Clone ID
    """
    orf81 = pd.read_csv(orf81_path)
    uorf = pd.read_csv(uorf_path)

    def _build(df, gene_col, id_col):
        ref = {}
        for _, row in df.iterrows():
            gene = str(row[gene_col]).strip()
            seq = str(row["nt_seq"]).strip()
            ref_id = row[id_col]
            ref.setdefault(gene, []).append((seq, ref_id))
        return ref

    orf81_ref = _build(orf81, "Gene_Symbol", "ORF_ID")
    uorf_ref = _build(uorf, "Gene Symbol", "Clone ID")

    log.info(f"ORFeome81 reference: {len(orf81_ref):,} genes, {len(orf81):,} entries")
    log.info(f"uORF reference:      {len(uorf_ref):,} genes, {len(uorf):,} entries")
    return orf81_ref, uorf_ref


# ---------------------------------------------------------------------------
# Library type detection
# ---------------------------------------------------------------------------


def library_for_barcode_len(barcode_len: int) -> str | None:
    """Returns 'orf81', 'uorf', or None if length is ambiguous."""
    if 28 <= barcode_len <= 32:
        return "orf81"
    if 22 <= barcode_len <= 26:
        return "uorf"
    return None


# ---------------------------------------------------------------------------
# Gene alias lookup via MyGene.info
# ---------------------------------------------------------------------------

_mg = mygene.MyGeneInfo()
_alias_cache: dict[str, list[str]] = {}


def get_aliases(gene_name: str) -> list[str]:
    """Returns a list of known aliases for gene_name via MyGene.info (cached)."""
    if gene_name in _alias_cache:
        return _alias_cache[gene_name]
    try:
        results = _mg.query(
            gene_name,
            fields="symbol,alias",
            species="human",
            size=1,
        )
        aliases = []
        if results.get("hits"):
            hit = results["hits"][0]
            if "alias" in hit:
                raw = hit["alias"]
                aliases = [raw] if isinstance(raw, str) else list(raw)
            if "symbol" in hit and hit["symbol"] != gene_name:
                aliases.insert(0, hit["symbol"])
        _alias_cache[gene_name] = aliases
        return aliases
    except Exception as e:
        log.debug(f"MyGene.info lookup failed for {gene_name}: {e}")
        _alias_cache[gene_name] = []
        return []


# ---------------------------------------------------------------------------
# Sequence matching
# ---------------------------------------------------------------------------


def find_best_match(
    query: str,
    candidates: list[tuple[str, object]],
    threshold: int,
) -> object | None:
    """Returns the ref_id of the closest candidate within threshold, or None."""
    best_dist = threshold + 1
    best_id = None
    for nt_seq, ref_id in candidates:
        dist = Levenshtein.distance(query, nt_seq, score_cutoff=threshold)
        if dist < best_dist:
            best_dist = dist
            best_id = ref_id
    return best_id


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------


def process():
    log.info(f"Input CSV:       {args.input}")
    log.info(f"Output CSV:      {args.output}")
    log.info(f"ORFeome81 ref:   {args.orf81_ref}")
    log.info(f"uORF ref:        {args.uorf_ref}")
    log.info(f"Threshold:       {args.threshold}")
    log.info(f"Log file:        {args.log_file}")

    orf81_ref, uorf_ref = load_references(args.orf81_ref, args.uorf_ref)

    df = pd.read_csv(args.input)
    log.info(f"Loaded {len(df):,} rows from {args.input}")

    df["barcode_len"] = df["barcode"].str.len()
    df["matched_ref_id"] = pd.NA

    stats = {
        "matched": 0,
        "no_library": 0,
        "gene_not_found": 0,
        "alias_resolved": 0,
        "no_seq_match": 0,
    }

    n_genes = df["gene_name"].nunique()
    for gene_name, group in tqdm.tqdm(
        df.groupby("gene_name"), total=n_genes, desc="Matching genes"
    ):
        # Determine library from most common barcode length in this group
        most_common_len = int(group["barcode_len"].mode().iloc[0])
        lib = library_for_barcode_len(most_common_len)

        if lib is None:
            stats["no_library"] += 1
            log.warning(
                f"Gene {gene_name}: barcode length {most_common_len} does not match "
                f"either library range (24±2 or 30±2) — skipping"
            )
            continue

        ref = orf81_ref if lib == "orf81" else uorf_ref

        # Look up gene; fall back to aliases if not found directly
        candidates = ref.get(gene_name)
        resolved_name = gene_name

        if candidates is None:
            aliases = get_aliases(gene_name)
            for alias in aliases:
                if alias in ref:
                    candidates = ref[alias]
                    resolved_name = alias
                    stats["alias_resolved"] += 1
                    log.debug(
                        f"Gene {gene_name}: not in {lib} reference, resolved via alias '{alias}'"
                    )
                    break

        if candidates is None:
            stats["gene_not_found"] += 1
            log.debug(
                f"Gene {gene_name}: not found in {lib} reference (including aliases)"
            )
            continue

        # Match each unique original_cdna against reference candidates once,
        # then map results back to all rows sharing that sequence
        seq_to_ref_id: dict[str, object] = {}
        for orig_seq in group["original_cdna"].dropna().unique():
            match = find_best_match(orig_seq, candidates, args.threshold)
            if match is not None:
                seq_to_ref_id[orig_seq] = match
                stats["matched"] += 1
            else:
                stats["no_seq_match"] += 1
                log.debug(
                    f"Gene {gene_name} ({resolved_name}): no reference match within "
                    f"threshold {args.threshold} for sequence of length {len(orig_seq)}"
                )

        if seq_to_ref_id:
            mask = group["original_cdna"].isin(seq_to_ref_id)
            idx = group.index[mask]
            df.loc[idx, "matched_ref_id"] = (
                group.loc[mask, "original_cdna"].map(seq_to_ref_id).values
            )

    df.drop(columns=["barcode_len"], inplace=True)
    df.to_csv(args.output, index=False)

    divider = "=" * 40
    thin = "-" * 40
    log.info(divider)
    log.info("REFERENCE ORF MATCHING SUMMARY")
    log.info(divider)
    log.info(f"Input CSV:                    {args.input}")
    log.info(f"Output CSV:                   {args.output}")
    log.info(f"Levenshtein threshold:        {args.threshold}")
    log.info(thin)
    log.info(f"Total rows:                   {len(df):,}")
    log.info(f"Unique sequences matched:     {stats['matched']:,}")
    log.info(f"Unique sequences unmatched:   {stats['no_seq_match']:,}")
    log.info(thin)
    log.info(f"Genes — no library match:     {stats['no_library']:,}")
    log.info(f"Genes — not in reference:     {stats['gene_not_found']:,}")
    log.info(f"Genes — resolved via alias:   {stats['alias_resolved']:,}")
    log.info(divider)


if __name__ == "__main__":
    process()
