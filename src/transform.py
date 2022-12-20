import tempfile
from collections import defaultdict
from os import makedirs
from os.path import abspath, basename, isfile, join
from re import finditer, sub
from typing import Dict, Iterable, List, Optional, Tuple
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.SeqIO import parse
from pileup_to_consensus import pileup_to_consensus
from reads_to_pileup import reads_to_pileup

from utils import (darken_color, file, get_sample_name_and_extenstion, pushd,
                   run, string_to_color, write_fasta)



# Extract the exons from the gene sequence. Take shortcuts because
# we have a very specific GFF file. The GFF and reference were
# transformed so everything is on plus strand, and seqid is the gene.
# So we don't have to worryabout phase or attributes.


def consensus_to_exons(consensus_path: str, gff: str, out_dir: str) -> str:

    sample_name, sample_ext = get_sample_name_and_extenstion(
        consensus_path, 'fasta')

    # Find the exon positions
    exon_positions = defaultdict(set)
    with file(gff) as f:
        for line in f:
            # Assumes everything's on the positive strand
            seqid, _, type_, start, end, _, strand, _, _ = line.strip().split('\t')
            assert strand == '+', 'Only + strand supported'
            if type_ == 'exon':
                exon_positions[seqid].update(
                    set(range(int(start) - 1, int(end))))

    consensus_exons_path = join(out_dir, f'{sample_name}_exons{sample_ext}')

    # Write the gene bases only in the exon positions
    consensus_exons = {
        r.id: ''.join(r.seq[i] for i in sorted(exon_positions[r.id]))
        for r in parse(consensus_path, 'fasta')
    }
    write_fasta(consensus_exons, consensus_exons_path, sort=True)

    # Let the caller know where to find the consensus exons
    return consensus_exons_path


def exons_to_exons_concat(exons_path: str, out_dir: str) -> str:
    sample_name, sample_ext = get_sample_name_and_extenstion(
        exons_path, 'fasta')
    exons_concat_path = join(out_dir, f'{sample_name}_concat{sample_ext}')
    if sample_name.endswith('_exons'):
        sample_name = sample_name[:-len('_exons')]
    all_genes_exons = ''.join(str(r.seq) for r in parse(exons_path, 'fasta'))
    write_fasta({sample_name: all_genes_exons}, exons_concat_path, sort=True)

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
    threads=1,
    trim=True,
    max_read_length=None,
) -> str:
    pileup_path = reads_to_pileup(
        fastq=fastq,
        reference=reference,
        out_dir=out_dir,
        threads=threads,
        trim=trim,
        max_read_length=max_read_length,
    )
    sample_name, _ = get_sample_name_and_extenstion(pileup_path, 'pileup')
    consensus_path = join(out_dir, f'{sample_name}.fasta')
    pileup_to_consensus(
        pileup_path=pileup_path,
        ref_path=reference,
        out_path=consensus_path,
        min_snp_depth=min_snp_depth,
        min_ref_depth=min_match_depth,
        hetero_min=hetero_min,
        hetero_max=hetero_max,
    )
    exons = consensus_to_exons(consensus_path, gff, out_dir)
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

def seq_to_pct_coverage(seq: str) -> float:
    seq = str(seq).upper()
    n_unknown = seq.count('?') + seq.count('N')
    pct_unknown = 100 * n_unknown / len(seq)
    return 100 - pct_unknown

def consensus_to_coverage(consensus, out_dir, step=1):
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(consensus, 'fasta')
    out_path = join(out_dir, f'{sample_name}.csv')
    cov_pcts = [seq_to_pct_coverage(r.seq) for r in parse(consensus, 'fasta')]

    pct_coverage_thresholds = list(range(0, 101, step))
    number_of_genes_with_pct_coverage_ge_thresholds = [
        len([pct for pct in cov_pcts if pct >= min_pct])
        for min_pct in pct_coverage_thresholds
    ]
    pd.DataFrame({
        'pct_exon_positions_covered': pct_coverage_thresholds,
        'number_of_genes': number_of_genes_with_pct_coverage_ge_thresholds,
    }).to_csv(out_path, index=None, header=None)
    return out_path

# This is a work around to ensure that multiqc can create the plot specified in
# config/multiqc_config.yaml.
# Unless there is at least one _mqc.yaml in the directory, multiqc will not 
# pick up any of the custom_data specified in the --config file. To ensure it
# does pick up those files, we create a file for each sample with only id.
def create_empty_config_required_for_gene_coverage_mqc(sample: str, out_dir: str):
    config_name = f'{sample}_empty_config_required_for_gene_coverage_mqc'
    with open(join(out_dir, f'{config_name}.yaml'), 'w') as f:
        f.write(f'id: "{config_name}"' + '\n')

def sample_report(sample_dir: str, sample_name: str):
    out_dir = join(sample_dir, 'report')
    for fastq_ext in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']:
        fastq_path = join(sample_dir, f'{sample_name}{fastq_ext}')
        if isfile(fastq_path):
            reads_to_fastqc(fastq_path, out_dir)
            break
    alignment_to_flagstat(join(sample_dir, f'{sample_name}.bam'), out_dir)
    consensus_to_coverage(join(sample_dir, f'{sample_name}.fasta'), out_dir)
    create_empty_config_required_for_gene_coverage_mqc(sample_name, out_dir)

def consensuses_to_coverage_table(
    consensus_paths: Iterable[str],
    collection_name: str,
    out_dir: str
) -> str:
    makedirs(out_dir, exist_ok=True)
    coverage_rows = []
    for consensus_path in consensus_paths:
        sample_name, _ = get_sample_name_and_extenstion(
            consensus_path, 'fasta')
        pct_coverages = []
        for r in parse(consensus_path, 'fasta'):
            pct_missing = (r.seq.count('?') + r.seq.count('N')) / len(r.seq)
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


def exons_concat_to_newick(exons_concat: str, out_dir: str, n_threads=1) -> str:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(exons_concat, 'fasta')
    # need absolute path because we're about to run raxml from the output directory
    exons_concat = abspath(exons_concat)

    # output name cannot contain slashes, so run raxml from the output directory
    with pushd(out_dir):
        # make a temporary file to prevent raxml log being written to screen
        with tempfile.TemporaryFile() as log:
            try:
                args = [
                    'raxmlHPC-PTHREADS-SSE3',
                     '-s', exons_concat,
                     '-m', 'GTRGAMMA',
                     '-n', f'{collection_name}.newick',
                     '-p', '100',
                ]
                if n_threads >= 2:
                    args += ['-T', str(n_threads)]
                run(args, out=log)
            except:
                # if raxml failed then raise an exception with the error message
                # which was written to the temporary file (redirected from stdout)
                log.seek(0)
                error = log.read().decode()
                # When run with a fasta file, RAxML reports it cannot be parsed as a phylip.
                # This message can make it harder to find the actual error if RAxML crashes,
                # so supress that message.
                error = error.replace(
                    "\nRAxML can't, parse the alignment file as phylip file \nit will now try to parse it as FASTA file\n\n",
                    ''
                )
                raise Exception(error)

    return join(out_dir, f'RAxML_bestTree.{collection_name}.newick')


def newick_to_imgs(newick_path: str, metadata_path: str, out_dir: str, img_fmt='pdf') -> Dict[str, str]:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(newick_path, 'newick')

    # Load metadata and styles
    tree = Phylo.read(newick_path, 'newick')
    tree.root_at_midpoint()
    tree.ladderize(reverse=True)
    sheets = pd.read_excel(metadata_path, sheet_name=None, engine='openpyxl')
    metadata = sheets['metadata'].fillna('?').astype(
        str).set_index('tree_name').to_dict()
    style = {}
    for sheet_name, sheet in sheets.items():
        if sheet_name == 'metadata':
            continue
        style[sheet_name] = sheet.set_index(
            sheet_name)[['color', 'marker']].to_dict(orient='index')

    show_labels = True
    n_leaves = len(tree.get_terminals())
    label_col = 'tree_new_name'

    img_out_paths = {}

    for style_col in list(style) + ['region']:

        img_out_path = join(out_dir, f'{collection_name}_{style_col}.{img_fmt}')
        img_out_paths[style_col] = img_out_path

        height = max(5, n_leaves / 6)
        width = 14
        fontsize = 11 - width / 5
        fig, ax = plt.subplots(1, 1, figsize=(width, height))
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

            default_style = {'color': string_to_color(
                style_val), 'marker': '●'}

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
                ax.text(x + 1.8 * font_size_x_units, y,
                        s, va='center', fontsize=fontsize)

        # Fill in the contiguous regions where the style_val is the same
        right = xmax + .15 * (xmax - xmin)
        text_right = xmax + .142 * (xmax - xmin)
        polygons = []
        to_fill = pd.DataFrame(name_to_pos).T.rename(columns={0: 'x', 1: 'y'})
        to_fill = to_fill.join(pd.DataFrame(name_to_style).T)
        to_fill['style_val'] = [metadata[style_col].get(
            s, '?') for s in to_fill.index]
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
        mean_contiguous_y = to_fill.groupby(
            ['style_val_index', 'style_val', 'color']).y.mean().to_frame().reset_index()
        for style_val, y, color in mean_contiguous_y[['style_val', 'y', 'color']].values:
            ax.text(text_right, y, style_val, color=darken_color(color),
                    ha='right', va='center', alpha=.5, fontsize=fontsize)

        # Make the branches thinner
        for collection in ax.collections:
            if list(collection.get_linewidths()) == [1.5]:
                collection.set_linewidths([0.5])

        # Add a legend
        for style_val in sorted(seen_style_vals):
            default_style = {'color': string_to_color(
                style_val), 'marker': '●'}
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
        ax.legend(title=' '.join(style_col.split('_')).title(), loc='upper left')

        # Hide the right and top spines
        for side in ['right', 'top']:
            ax.spines[side].set_visible(False)

        # Make the left and top axes less prominent
        faint_color = '#BBB'
        for side in ['left', 'bottom']:
            ax.spines[side].set_color(faint_color)
        ax.tick_params(axis='x', colors=faint_color)
        ax.tick_params(axis='y', colors=faint_color)
        ax.yaxis.label.set_color(faint_color)
        ax.xaxis.label.set_color(faint_color)

        ax.set_ylabel(None)

        ax.set_xlim(xmin, right)

        plt.tight_layout()
        plt.savefig(img_out_path)
        plt.close()

    return img_out_paths

def reads_list_to_exons_concat_with_report(
    fastq_paths: List[str],
    reference: str,
    gff: str,
    out_dirs: List[str],
    multiqc_config: str,
    threads=1,
    trim=True,
):
    for fastq_index, (fastq, out_dir) in enumerate(zip(fastq_paths, out_dirs)):
        sample_name = get_sample_name_and_extenstion(fastq, 'fastq')[0]
        print(f'{fastq_index + 1}/{len(fastq_paths)} {sample_name}:', end=' ', flush=True)
        reads_to_exons_concat(
            fastq=fastq,
            reference=reference,
            gff=gff,
            out_dir=out_dir,
            threads=threads,
            trim=trim,
        )
        print('assessing', flush=True)
        sample_report(out_dir, sample_name)
    print('Report: compiling')
    report_dirs = [join(out_dir, 'report') for out_dir in out_dirs]
    run(['multiqc', '--config', multiqc_config] + report_dirs, out='/dev/null')


def exon_concat_paths_to_tree_input(
    exon_concat_paths: List[str],
    starting_tree_input: str,
    new_tree_input: str,
) -> str:
    with open(new_tree_input, 'w') as f_out:
        for path in [*exon_concat_paths, starting_tree_input]:
            with file(path) as f_in:
                for line in f_in:
                    f_out.write(line)
    # You already know where to find it, but for consistency
    return new_tree_input


def exon_concat_paths_to_tree_imgs(
    exon_concat_paths: List[str],
    starting_tree_input: str,
    tree_name: str,
    out_dir: str,
    metadata_path: str,
    n_threads=1,
    img_fmt='pdf',
) -> Dict[str, str]:
    makedirs(out_dir, exist_ok=True)
    tree_input = join(out_dir, f'{tree_name}.fasta')
    exon_concat_paths_to_tree_input(
        exon_concat_paths=exon_concat_paths,
        starting_tree_input=starting_tree_input,
        new_tree_input=tree_input,
    )
    print('Tree: making', end=' ', flush=True)
    newick = exons_concat_to_newick(
        exons_concat=tree_input,
        out_dir=out_dir,
        n_threads=n_threads,
    )
    print('visualising', flush=True)
    imgs = newick_to_imgs(
        newick_path=newick,
        metadata_path=metadata_path,
        out_dir=out_dir,
        img_fmt=img_fmt,
    )
    return imgs
