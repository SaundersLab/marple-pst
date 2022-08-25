
from Bio.SeqIO import parse
import pandas as pd
from utils import file, get_sample_name_and_extenstion, run, pushd
from os.path import abspath, join
from os import makedirs
from re import finditer, sub
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import tempfile

# Make a pileup and return the path to it
def reads_to_pileup(fastq: str, reference: str, out_dir: str) -> str:
    makedirs(out_dir, exist_ok=True)

    sample_name, sample_ext = get_sample_name_and_extenstion(fastq, 'fastq')

    porechopped = join(out_dir, f'{sample_name}_porechopped{sample_ext}')
    run(['porechop', '-i', fastq], porechopped)

    run(['bwa', 'index', reference])

    aligned_sam = join(out_dir, f'{sample_name}.sam')
    run(['bwa', 'mem', reference, porechopped], aligned_sam)

    aligned_bam = join(out_dir, f'{sample_name}.bam')
    run(['samtools', 'view', '-S', '-b', aligned_sam], aligned_bam)

    sorted_bam = join(out_dir, f'{sample_name}_sorted.bam')
    run(['samtools', 'sort', aligned_bam], sorted_bam)

    run(['samtools', 'faidx', reference])

    pileup = join(out_dir, f'{sample_name}.pileup')
    run(['samtools', 'mpileup', '-f', reference, sorted_bam], pileup)

    # Let the caller know where to find the pileup file
    return pileup

# Make a consensus and return the path to it
def pileup_to_consensus(
    pileup_path: str,
    reference_path: str,
    out_dir: str,
    min_snp_depth: int,
    min_match_depth: int,
    hetero_min=.25,
    hetero_max=.75,
) -> str:

    sample_name, sample_ext = get_sample_name_and_extenstion(pileup_path, 'pileup')

    snp_ratios_path = join(out_dir, f'{sample_name}_snp_ratios.tsv')
    pileup_to_snp_ratios(pileup_path, snp_ratios_path)

    snp_freq_path = join(out_dir, f'{sample_name}_snp_freq.tsv')
    snp_ratios_to_snp_freq(snp_ratios_path, snp_freq_path, min_snp_depth, hetero_min, hetero_max)

    snp_freq = pd.read_csv(snp_freq_path, sep='\t', header=None)
    snp_freq.columns = ['gene', 'pos', 'pos1', 'ref', 'depth', 'allele', 'genotype', 'ratios']
    # Drop the position where 3+ bases qualified (e.g. T ACT ? 0.2,0.6,0.2)
    snp_freq = snp_freq[snp_freq.genotype != '?']

    snp_ratios = pd.read_csv(snp_ratios_path, sep='\t', header=None)
    snp_ratios.columns = ['gene', 'pos', 'ref', 'depth', 'alts', 'ratios']

    allele_to_code = {
        'AA': 'A',  'CC': 'C',  'GG': 'G',  'TT': 'T',  'UU': 'U',
        'AT': 'W',  'CG': 'S',  'AC': 'M',  'GT': 'K',  'AG': 'R',  'CT': 'Y',
        'TA': 'W',  'GC': 'S',  'CA': 'M',  'TG': 'K',  'GA': 'R',  'TC': 'Y',
    }

    # Start with an empty consensus
    consensus = {r.id: ['N'] * len(r.seq) for r in parse(reference_path, 'fasta')}

    # Fill in the consensus using alternate bases at positions with coverage >= min_depth
    for gene, pos, genotype, allele in snp_freq[['gene', 'pos', 'genotype', 'allele']].values:
        consensus[gene][pos] = allele_to_code.get(allele, consensus[gene][pos])

    # Fill in the consensus using matches to reference with coverage >= min_match_depth
    ref_matches = snp_ratios[(snp_ratios.alts == snp_ratios.ref) & (snp_ratios.depth >= min_match_depth)]
    for gene, pos, ref in ref_matches[['gene', 'pos', 'ref']].values:
        consensus[gene][pos] = ref

    # Write out the consensus
    consensus_path = join(out_dir, f'{sample_name}.fasta')
    with file(consensus_path, 'wt') as f_out:
        for gene in sorted(consensus):
            f_out.write('>' + gene + '\n' + ''.join(consensus[gene]) + '\n')

    # Let the caller know where to find the consensus file
    return consensus_path

def _subtract_indels_from_counts(upper_reads: str, base_counts: Dict[str, int]) -> None:
    for indel in finditer(r'([+|-])(\d+)(\w+)', upper_reads):
        indel_size = int(indel.group(2))
        indel_sequence = indel.group(3)[:indel_size]
        for base in base_counts:
            base_counts[base] -= indel_sequence.count(base)

def _base_ratios_from_reads(reads: str, depth: int, ref: str) -> Dict[str, float]:
    upper_reads = reads.upper()
    
    # ^G. marks the start of a read segment where the ASCII minus 33 of G is the mapping quality
    # and the character AFTER those 2 is the actual read - in this case a match to the reference. 
    # Drop the first 2 characters as they are not actual reads
    if '^' in upper_reads:
        upper_reads = sub('\^.', '', upper_reads)
    
    base_counts = {base: upper_reads.count(base) for base in 'ACTG'}
    if ref in 'ACTG':
        base_counts[ref] += upper_reads.count('.') + upper_reads.count(',')
    if '+' in upper_reads or '-' in upper_reads:
        _subtract_indels_from_counts(upper_reads, base_counts)
    # Only calculate ratios for present bases
    base_ratios = {
        base: round(base_counts[base] / depth, 3)
        for base in base_counts if base_counts[base]
    }
    return base_ratios

def _pileup_row_to_snp_ratio_row(row: str) -> str:
    contig, pos_str, ref, depth_str, reads, *_ = row.split('\t')
    depth = int(depth_str)
    base_to_ratio = _base_ratios_from_reads(reads, depth, ref)
    # sort by base for backward compatability
    base_ratios = sorted(base_to_ratio.items(), key=lambda base_and_ratio: base_and_ratio[0])
    bases_str = ','.join(base for base, ratio in base_ratios)
    ratios_str = ','.join(str(ratio) for base, ratio in base_ratios)
    return '\t'.join((contig, pos_str, ref, depth_str, bases_str, ratios_str))

def pileup_to_snp_ratios(pileup_path: str, snp_ratios_path: str) -> None:
    with file(pileup_path, 'rt') as f_in, file(snp_ratios_path, 'wt') as f_out:
        for row in f_in:
            f_out.write(_pileup_row_to_snp_ratio_row(row) + '\n')

# Determine the genotype at each position and filter by read coverage
def _get_genotype_and_valid_bases_and_valid_ratios(
    ref: str, bases: List[str], ratio_strs: list, hetero_min: float, hetero_max: float
) -> Optional[Tuple[str, str, str]]:

    ratios = map(float, ratio_strs)

    genotype = ''
    valid_bases = ''
    valid_ratios = ''

    for (base, ratio, ratio_str) in zip(bases, ratios, ratio_strs):
        if ratio < hetero_min:
            continue
        valid_ratios += ratio_str + ','
        if ratio >= hetero_max:
            genotype = '0/0' if base == ref else '1/1'
            valid_bases = base + base
            break
        valid_bases += base

    if not valid_bases:
        return None

    if not genotype:
        if len(valid_bases) == 1:
            if valid_bases == ref:  # If only ref has min <= freq <= max
                return None
            genotype = '1/1'
            valid_bases = valid_bases + valid_bases
        elif len(valid_bases) == 2:
            genotype = '0/1' if ref in valid_bases else '1/2'
        else:
            genotype = '?'

    valid_ratios = valid_ratios.rstrip(',')
    return genotype, valid_bases, valid_ratios

def _snp_ratio_row_to_snp_freq_row(
    row: str, min_depth: int, hetero_min: float, hetero_max: float
) -> Optional[str]:

    contig, pos, ref, depth, bases_str, ratios = row.split('\t')

    if not bases_str or int(depth) < min_depth:
        return None

    bases = bases_str.split(',')
    ratio_strs = ratios.rstrip().split(',')

    genotype_and_valid_bases_and_valid_ratios = _get_genotype_and_valid_bases_and_valid_ratios(
        ref, bases, ratio_strs, hetero_min, hetero_max
    )

    if not genotype_and_valid_bases_and_valid_ratios:
        return None

    genotype, valid_bases, valid_ratios = genotype_and_valid_bases_and_valid_ratios

    return '\t'.join((contig, str(int(pos) - 1), pos, ref, depth, valid_bases, genotype, valid_ratios))

def snp_ratios_to_snp_freq(
    snp_ratios_path: str,
    snp_freq_path: str,
    min_depth: int,
    hetero_min: float,
    hetero_max: float,
    ) -> None:
    with file(snp_ratios_path, 'rt') as f_in, file(snp_freq_path, 'wt') as f_out:
        for row in f_in:
            snp_freq_row = _snp_ratio_row_to_snp_freq_row(
                row, min_depth, hetero_min, hetero_max
            )
            if snp_freq_row:
                f_out.write(snp_freq_row + '\n')

# Extract the exons from the gene sequence. Take shortcuts because
# we have a very specific GFF file. The GFF and reference were 
# transformed so everything is on plus strand, and seqid is the gene.
# So we don't have to worryabout phase or attributes.
def consensus_to_exons(consensus_path: str, gff: str, out_dir: str) -> str:

    sample_name, sample_ext = get_sample_name_and_extenstion(consensus_path, 'fasta')

    # Find the exon positions
    exon_positions = defaultdict(set)
    for line in file(gff):
        # Assumes everything's on the positive strand
        seqid, source, type_, start, end, score, strand, phase, attributes = line.strip().split('\t')
        assert strand == '+', 'Only + strand supported'
        if type_ == 'exon':
            exon_positions[seqid].update(set(range(int(start) - 1, int(end) - 1)))

    consensus_exons_path = join(out_dir, f'{sample_name}_exons{sample_ext}')

    # Write the gene bases only in the exon positions
    with file(consensus_exons_path, 'wt') as cds_out:
        for r in parse(consensus_path, 'fasta'):
            exons_seq = ''.join(r.seq[i] for i in sorted(exon_positions[r.id]))
            cds_out.write('>' + r.id + '\n' + exons_seq + '\n')

    # Let the caller know where to find the consensus exons
    return consensus_exons_path

def exons_to_exons_concat(exons_path: str, out_dir: str) -> str:
    sample_name, sample_ext = get_sample_name_and_extenstion(exons_path, 'fasta')
    exons_concat_path = join(out_dir, f'{sample_name}_concat{sample_ext}')
    all_gene_exon_seqs = ''.join(str(r.seq) for r in parse(exons_path, 'fasta'))
    with file(exons_concat_path, 'wt') as f_out:
        f_out.write('>' + sample_name + '\n' + all_gene_exon_seqs + '\n')
    return exons_concat_path

def reads_to_exons_concat(
    fastq: str,
    reference: str,
    gff: str,
    out_dir: str,
    min_snp_depth: int = 20,
    min_match_depth: int = 2,
    hetero_min: float = .2,
    hetero_max: float = .8,
) -> str:

    pileup = reads_to_pileup(fastq, reference, out_dir)
    consensus = pileup_to_consensus(
        pileup, reference, out_dir, min_snp_depth, min_match_depth, hetero_min, hetero_max
    )
    exons = consensus_to_exons(consensus, gff, out_dir)
    exons_concat = exons_to_exons_concat(exons, out_dir)
    return exons_concat

def reads_to_fastqc(fastq: str, out_dir: str) -> Tuple[str, str]:
    makedirs(out_dir, exist_ok=True)
    sample_name, sample_ext = get_sample_name_and_extenstion(fastq, 'fastq')
    fastqc_page = join(out_dir, f'{sample_name}_fastqc.html')
    fastqc_data = join(out_dir, f'{sample_name}_fastqc.zip')
    run(['fastqc', '--quiet', '-o', out_dir, fastq])
    return fastqc_page, fastqc_data

def alignment_to_flagstat(alignment: str, out_dir: str) -> str:
    makedirs(out_dir, exist_ok=True)
    sample_name, sample_ext = get_sample_name_and_extenstion(alignment, 'alignment')
    flagstat = join(out_dir, f'{sample_name}.txt')
    run(['samtools', 'flagstat', alignment], flagstat)
    return flagstat

def exons_concat_to_newick(exons_concat: str, out_dir: str, n_threads=1) -> str:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(exons_concat, 'fasta')
    # need absolute path because we're about to run raxml from the output directory
    exons_concat = abspath(exons_concat)
    
    # output name cannot contain slashes, so run raxml from the output directory
    with pushd(out_dir):
        # make a temporary file to prevent raxml log being written to screen
        with tempfile.TemporaryFile(mode='wt') as log:
            try:
                run(['raxmlHPC-PTHREADS-SSE3',
                    '-T', str(n_threads),
                    '-s', exons_concat,
                    '-m', 'GTRGAMMA',
                    '-n', collection_name,
                    '-p', '100',
                ], out=log)
            except:
                # if raxml failed then raise an exception with the error message
                # which was written to the temporary file (redirected from stdout)
                log.seek(0)
                error = log.read()
                # supress warning that looks like an error
                error = error.replace("\nRAxML can't, parse the alignment file as phylip file \nit will now try to parse it as FASTA file\n\n", '')
                raise Exception(error)

    return join(out_dir, f'RAxML_bestTree.{collection_name}.newick')

