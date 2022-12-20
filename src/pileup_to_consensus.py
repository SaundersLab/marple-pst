from re import sub, search
from typing import Tuple, Dict, Union, List
from Bio.SeqIO import parse
import argparse
from utils import file

def base_ratios_from_reads(
    ref: str, depth: int, reads: str,
) -> Dict[str, float]:

    # Don't count mapping qualities at the start of read segments as bases
    # e.g. ^G. is a read segment start and the ASCII minus 33 of G is the
    # mapping quality. Drop the 2 characters that are not reads.
    if '^' in reads:
        reads = sub('\^.', '', reads)

    # Remove pileup indels so you don't count the bases in an indel e.g. +3ATT
    if '+' in reads or '-' in reads:
        while match := search(r'[+-](\d+)', reads):
            reads = reads[:match.start()] + reads[match.end() + int(match[1]):]

    upper_reads = reads.upper()
    base_counts = {base: upper_reads.count(base) for base in 'ACTG'}

    upper_ref = ref.upper()
    if upper_ref in 'ACTG':
        base_counts[upper_ref] += reads.count('.') + reads.count(',')

    # Only calculate ratios for present bases
    return {base: n / depth for base, n in base_counts.items() if n}


def get_genotype_and_valid_base_ratios(
    ref: str,
    base_ratios: Dict[str, float],
    hetero_min: float,
    hetero_max: float,
) -> Union[Tuple[str, Dict[str, float]], Tuple[None, None]]:
    """
    Use the heterozygosity thresholds to filter the base ratios and determine the genotype:
    Genotype  Condition	                                        Example
    0/0	      pref ≥ hetero_max or only refbase qualified	    G GG 0/0 1
    0/1	      ref base and 1 alt base qualified	                T CT 0/1 0.217,0.783
    1/1	      palt 1 ≥ hetero_max or only 1 alt base qualified	T CC 1/1 0.822
    1/2	      2 alt bases qualified and ref base did not	    C AT 1/2 0.476,0.524
    ?	      3+ bases qualified	                            T ACT ? 0.2,0.6,0.2
    None      No bases qualify                                    
    """
    assert hetero_max > hetero_min

    valid_base_ratios = {}
    for base, ratio in base_ratios.items():
        if ratio >= hetero_min:
            valid_base_ratios[base] = ratio
            if ratio >= hetero_max:
                return '0/0' if base == ref else '1/1', valid_base_ratios

    if not valid_base_ratios:
        return None, None

    # If only ref has min <= freq <= max
    if list(valid_base_ratios) == [ref]:
        return None, None

    if len(valid_base_ratios) == 1:
        genotype = '1/1'
    elif len(valid_base_ratios) == 2:
        genotype = '0/1' if ref in valid_base_ratios else '1/2'
    else:
        genotype = '?'

    return genotype, valid_base_ratios


def pileup_row_to_consensus_at_locus(
    row: str,
    min_ref_depth: int = 2,
    min_snp_depth: int = 10,
    hetero_min: float = .2,
    hetero_max: float = .8,
) -> Tuple[str, int, str]:
    alleles_to_code = {
        'A': 'A',  'C': 'C',  'G': 'G',  'T': 'T',
        'AT': 'W',  'CG': 'S',  'AC': 'M',  'GT': 'K',  'AG': 'R',  'CT': 'Y',
        'TA': 'W',  'GC': 'S',  'CA': 'M',  'TG': 'K',  'GA': 'R',  'TC': 'Y',
    }
    contig, pos_str, ref, depth_str, reads, *_ = row.split('\t')
    depth = int(depth_str)
    pos = int(pos_str)
    null_consensus = contig, pos - 1, 'N'
    if depth < min_ref_depth:
        return null_consensus
    base_ratios = base_ratios_from_reads(ref, depth, reads)
    if list(base_ratios) == [ref]:
        return contig, pos - 1, ref
    if depth < min_snp_depth:
        return null_consensus
    genotype, valid_base_ratios = get_genotype_and_valid_base_ratios(
        ref, base_ratios, hetero_min, hetero_max
    )
    if valid_base_ratios is None or genotype == '?':
        return null_consensus
    alleles = ''.join(valid_base_ratios)
    return contig, pos - 1, alleles_to_code[alleles]


def pileup_to_consensus(
    pileup_path: str, ref_path: str, out_path: str, **kwargs,
) -> None:

    # Start with an empty consensus
    consensus = {r.id: ['N'] * len(r.seq) for r in parse(ref_path, 'fasta')}
    with file(pileup_path) as f:
        for row in f:
            contig, pos, consensus[contig][pos] = pileup_row_to_consensus_at_locus(
                row, **kwargs
            )

    with file(out_path, 'wt') as f:
        for contig in consensus:
            f.write('>' + contig + '\n' + ''.join(consensus[contig]) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a consensus genome from a pileup file.')
    parser.add_argument('--pileup', dest='pileup_path', type=str, required=True)
    parser.add_argument('--ref', dest='ref_path', help='Reference FASTA', type=str, required=True)
    parser.add_argument('--out', dest='out_path', help='Output consensus FASTA', type=str, required=True)
    parser.add_argument('--min_ref_depth', type=int, default=2)
    parser.add_argument('--min_snp_depth', type=int, default=10)
    parser.add_argument('--hetero_min', type=float, default=.2)
    parser.add_argument('--hetero_max', type=float, default=.8)
    args = vars(parser.parse_args())
    pileup_to_consensus(**args)
