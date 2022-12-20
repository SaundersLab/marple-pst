import tempfile
from os import makedirs
from os.path import abspath, isfile, join
from typing import Dict, Iterable, List, Tuple
import pandas as pd

from Bio.SeqIO import parse
from pileup_to_consensus import pileup_to_consensus
from reads_to_pileup import reads_to_pileup
from newick_to_images import newick_to_images

from report import report

from utils import (file, get_sample_name_and_extenstion, pushd,
                   run, write_fasta)

from consensus_to_exons import consensus_to_exons

def concat_fasta_sequences(
    fasta_path: str, concat_fasta_path: str, header: str
) -> None:
    """
    Concatenate all the sequences in a fasta file into a single record.
    """
    sequences_combined = ''.join(str(r.seq) for r in parse(fasta_path, 'fasta'))
    write_fasta({header: sequences_combined}, concat_fasta_path)

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
    exons_path = join(out_dir, f'{sample_name}_exons.fasta')
    consensus_to_exons(consensus_path, gff, exons_path)
    exons_concat_path = join(out_dir, f'{sample_name}_exons_concat.fasta')
    concat_fasta_sequences(exons_path, exons_concat_path, sample_name)
    return exons_concat_path



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
        report(out_dir, sample_name)
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
