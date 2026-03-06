import argparse
import pandas as pd
import tqdm
import pysam
from itertools import groupby
from collections import Counter

# Command line parser
# ----------------------------
parser = argparse.ArgumentParser(description="Generate ORF consensus sequences")
parser.add_argument(
    "-i",
    "--input",
    help="Input BAM file"
)
parser.add_argument(
    "-o",
    "--output",
    help="Output CSV file",
    default="raw_consensus_sequences.csv",
)
parser.add_argument(
    "-m",
    "--min-reads",
    type=int,
    help="Minimum reads per barcode to generate a consensus",
    default=3,
)
parser.add_argument(
    "-f",
    "--ref-fasta",
    help="Gencode Reference FASTA file",
)


# Configuration
# ----------------------------
args = parser.parse_args()
BAM_FILE = args.input
REF_FASTA = args.ref_fasta
MIN_READS = args.min_reads
OUTPUT_CSV = args.output


def get_bc_tag(read):
    return read.get_tag("BC") if read.has_tag("BC") else None

def get_pileup_consensus(reads, ref_len):
    """Generates consensus sequence using positional pileup."""
    pileup = [Counter() for _ in range(ref_len)]
    for read in reads:
        for q_pos, r_pos in read.get_aligned_pairs():
            if r_pos is not None and 0 <= r_pos < ref_len:
                base = read.query_sequence[q_pos] if q_pos is not None else '-'
                pileup[r_pos][base] += 1
    
    res = []
    for counts in pileup:
        if not counts: continue
        best = counts.most_common(1)[0][0]
        if best != '-': res.append(best)
    return "".join(res)

def process_to_csv():
    # Load reference FASTA for the "original" sequences
    print("Loading reference FASTA...")
    ref_seqs = pysam.FastaFile(REF_FASTA)
    
    results = []
    
    with pysam.AlignmentFile(BAM_FILE, "rb", check_sq=False) as bam:
        reader = bam.fetch(until_eof=True)
        
        print("Processing barcodes...")
        for bc, group in tqdm.tqdm(groupby(reader, key=get_bc_tag)):
            if not bc: continue
            
            reads = list(group)
            if len(reads) < MIN_READS: continue
            
            # Determine majority reference
            ref_counts = Counter(r.reference_name for r in reads if r.reference_name)
            if not ref_counts: continue
            full_header = ref_counts.most_common(1)[0][0]
            
            # 1. Barcode sequence (bc)
            # 2. Full GenCode header (full_header)
            # 3. GenCode cDNA sequence (from FASTA)
            try:
                original_cdna = ref_seqs.fetch(full_header)
            except KeyError:
                original_cdna = "REF_NOT_FOUND"
                
            # Sequences cDNA consensus sequence
            ref_len = reads[0].header.get_reference_length(full_header)
            consensus = get_pileup_consensus([r for r in reads if r.reference_name == full_header], ref_len)
            
            results.append({
                "barcode": bc,
                "gencode_header": full_header,
                "original_cdna": original_cdna,
                "consensus_cdna": consensus
            })

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Done! Results saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    process_to_csv()