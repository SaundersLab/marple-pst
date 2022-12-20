import tempfile
from collections import defaultdict
from os import makedirs
from os.path import abspath, basename, isfile, join
from re import finditer, sub
from typing import Dict, Iterable, List, Optional, Tuple
import subprocess
import pandas as pd

from Bio.SeqIO import parse
from pileup_to_consensus import pileup_to_consensus
from reads_to_pileup import reads_to_pileup
from newick_to_images import newick_to_images

from utils import (file, get_sample_name_and_extenstion, pushd,
                   run, write_fasta)



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




def reads_list_to_exons_concat_with_report(
    fastq_paths: List[str],
    reference: str,
    gff: str,
    out_dirs: List[str],
    multiqc_config: str,
    threads=1,
    trim=True,
    max_read_length=None,
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
            max_read_length=max_read_length,
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
    imgs = newick_to_images(
        newick_path=newick,
        metadata_path=metadata_path,
        out_dir=out_dir,
        img_fmt=img_fmt,
    )
    return imgs
