import tempfile
from collections import defaultdict
from os import makedirs
from os.path import abspath, join
from re import finditer, sub
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.SeqIO import parse

from utils import (darken_color, file, get_sample_name_and_extenstion, pushd,
                   run, string_to_color)


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

    sample_name, _ = get_sample_name_and_extenstion(pileup_path, 'pileup')

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
    for gene, pos, _, allele in snp_freq[['gene', 'pos', 'genotype', 'allele']].values:
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
    bases_str = ','.join(base for base, _ in base_ratios)
    ratios_str = ','.join(str(ratio) for _, ratio in base_ratios)
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
        seqid, _, type_, start, end, _, strand, _, _ = line.strip().split('\t')
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
    hetero_min: float = .25,
    hetero_max: float = .75,
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
    sample_name, _ = get_sample_name_and_extenstion(fastq, 'fastq')
    fastqc_page = join(out_dir, f'{sample_name}_fastqc.html')
    fastqc_data = join(out_dir, f'{sample_name}_fastqc.zip')
    run(['fastqc', '--quiet', '-o', out_dir, fastq])
    return fastqc_page, fastqc_data

def alignment_to_flagstat(alignment: str, out_dir: str) -> str:
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(alignment, 'alignment')
    flagstat = join(out_dir, f'{sample_name}.txt')
    run(['samtools', 'flagstat', alignment], flagstat)
    return flagstat

def consensuses_to_coverage_table(consensus_paths: Iterable[str], collection_name: str, out_dir: str) -> str:
    makedirs(out_dir, exist_ok=True)
    coverage_rows = []
    for consensus_path in consensus_paths:
        sample_name, _ = get_sample_name_and_extenstion(consensus_path, 'fasta')
        pct_coverages = []
        for r in parse(consensus_path, 'fasta'):
            length = len(r.seq)
            pct_missing = (r.seq.count('?') + r.seq.count('N')) / length
            pct_coverages.append(1 - pct_missing)
        coverage_rows.append({
            'sample': sample_name,
            '<50% coverage': len([p for p in pct_coverages if p < .5]),
            '≥50% coverage': len([p for p in pct_coverages if p >= .5]),
        })
    coverage_df = pd.DataFrame(coverage_rows)
    coverage_table_path = join(out_dir, f'{collection_name}.tsv')
    coverage_df.to_csv(coverage_table_path, sep='\t', index=None)
    return coverage_table_path

def exons_concat_to_newick(exons_concat: str, out_dir: str, n_threads=2) -> str:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(exons_concat, 'fasta')
    # need absolute path because we're about to run raxml from the output directory
    exons_concat = abspath(exons_concat)
    
    # output name cannot contain slashes, so run raxml from the output directory
    with pushd(out_dir):
        # make a temporary file to prevent raxml log being written to screen
        with tempfile.TemporaryFile() as log:
            try:
                run(['raxmlHPC-PTHREADS-SSE3',
                    '-T', str(n_threads),
                    '-s', exons_concat,
                    '-m', 'GTRGAMMA',
                    '-n', f'{collection_name}.newick',
                    '-p', '100',
                ], out=log)
            except:
                # if raxml failed then raise an exception with the error message
                # which was written to the temporary file (redirected from stdout)
                log.seek(0)
                error = log.read().decode()
                # supress warning that looks like an error
                error = error.replace("\nRAxML can't, parse the alignment file as phylip file \nit will now try to parse it as FASTA file\n\n", '')
                raise Exception(error)

    return join(out_dir, f'RAxML_bestTree.{collection_name}.newick')


def newick_to_pdfs(newick_path: str, metadata_path: str, out_dir: str) -> Dict[str, str]:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(newick_path, 'newick')

    # Load metadata and styles
    tree = Phylo.read(newick_path, 'newick')
    tree.ladderize()
    sheets = pd.read_excel(metadata_path, sheet_name=None, engine='openpyxl')
    metadata = sheets['metadata'].fillna('?').astype(str).set_index('tree_name').to_dict()
    style = {}
    for sheet_name, sheet in sheets.items():
        if sheet_name == 'metadata':
            continue
        style[sheet_name] = sheet.set_index(sheet_name)[['color', 'marker']].to_dict(orient='index')

    show_labels = True
    n_leaves = len(tree.get_terminals())
    label_col = 'tree_new_name'

    pdf_out_paths = {}

    for style_col in list(style) + ['region']:

        pdf_out_path = join(out_dir, f'{collection_name}_{style_col}.pdf')
        pdf_out_paths[style_col] = pdf_out_path

        size = max(12, n_leaves / 12)
        fontsize = 11 - size / 5
        fig, ax = plt.subplots(figsize=(.9 * size, size))
        ymin, ymax = (0, n_leaves)
        ax.set_ylim(ymin, ymax)

        Phylo.draw(tree, axes=ax, do_show=False)

        xmin, xmax = ax.get_xlim()

        # Get dimensions of axis in pixels
        x1, x2 = ax.get_window_extent().get_points()[:, 0]
        # Get unit scale
        xscale = (x2 - x1) / (xmax-xmin)
        # Get width of font in data units
        font_size_x_units = fontsize / xscale 

        seen_style_vals = set()
        
        name_to_pos = {}
        name_to_style = {}
        
        # Update the markers and labels
        texts = [t for t in ax.texts]
        for t in texts:
            s = t.get_text().strip()
            
            style_val = metadata[style_col].get(s, '?')
            seen_style_vals.add(style_val)

            default_style = {'color': string_to_color(style_val), 'marker': '●'}

            if (style_val == '?') or (style_col not in style):
                val_style = default_style
            else:
                val_style = style[style_col].get(style_val, default_style)
                
            name_to_pos[s] = t.get_position()
            name_to_style[s] = val_style

            color = val_style['color']
            marker = val_style['marker']

            t.set_text(marker)
            t.set_color(color)

            t.set_size(fontsize * 1.2)
            x, y = t.get_position()
            if show_labels:
                s = metadata[label_col].get(s, s)
                ax.text(x + 1.2 * font_size_x_units, y, s, va='center', fontsize=fontsize) 
            
        # Fill in the contiguous regions where the style_val is the same
        right = xmax
        polygons = []
        to_fill = pd.DataFrame(name_to_pos).T.rename(columns={0: 'x', 1: 'y'})
        to_fill = to_fill.join(pd.DataFrame(name_to_style).T)
        to_fill['style_val'] = [metadata[style_col].get(s, '?') for s in to_fill.index]
        to_fill = to_fill.sort_values(by='y')
        poly_x = []
        poly_y: list = []    
        curr_style_val = list(to_fill.style_val)[0]
        curr_color = list(to_fill.color)[0]
        for x, y, color, style_val in to_fill[['x', 'y', 'color', 'style_val']].values:
            if style_val != curr_style_val and poly_y:
                poly_x += [right, right]
                poly_y += [poly_y[-1], poly_y[0]]
                polygons += [poly_x, poly_y, curr_color]
                poly_x = []
                poly_y = []
                curr_style_val = style_val
                curr_color = color
            poly_x += [x + 1.2 * font_size_x_units] * 2
            poly_y += [y - .5, y + .5]
        poly_x += [right, right]
        poly_y += [y + .5, poly_y[0]]
        polygons += [poly_x, poly_y, curr_color]
        ax.fill(*polygons, alpha=.15, ec='#555', lw=.25)
        
        # Add a label to each filled region
        style_val_index = 0
        curr_style_val = None
        style_val_indices = []
        for style_val in to_fill.style_val:
            if style_val != curr_style_val:
                style_val_index += 1
                curr_style_val = style_val
            style_val_indices.append(style_val_index)
        to_fill['style_val_index'] = style_val_indices
        mean_contiguous_y = to_fill.groupby(['style_val_index', 'style_val', 'color']).y.mean().to_frame().reset_index()
        for style_val, y, color in mean_contiguous_y[['style_val', 'y', 'color']].values:
            ax.text(right, y, style_val, color=darken_color(color), ha='right', va='center', alpha=.5, fontsize=fontsize)
          
        # Make the branches thinner
        for collection in ax.collections:
            if list(collection.get_linewidths()) == [1.5]:
                collection.set_linewidths([0.5])

        # Add a legend
        for style_val in sorted(seen_style_vals):
            default_style = {'color': string_to_color(style_val), 'marker': '●'}
            if (style_val == '?') or (style_col not in style):
                val_style = default_style
            else:
                val_style = style[style_col].get(style_val, default_style)
            color = val_style['color']
            scatter_style = {
                '●': {'marker': 'o', 'fc': color, 'ec': color},
                '■': {'marker': 's', 'fc': color, 'ec': color},
                '○': {'marker': 'o', 'fc': '#FFF', 'ec': color},
            }[val_style['marker']]
            ax.scatter([], [], label=style_val, **scatter_style)
        ax.legend(title=' '.join(style_col.split('_')).title(), loc=(1.04,0))

        # Hide the right and top spines
        for side in ['right', 'top']:
            ax.spines[side].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(pdf_out_path)
        plt.close()

    return pdf_out_paths
