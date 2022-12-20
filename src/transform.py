from os import makedirs
from os.path import join
from typing import Dict, Iterable, List
import pandas as pd

from Bio.SeqIO import parse
from newick_to_images import newick_to_images
from fasta_to_newick import fasta_to_newick
from report import report

from utils import file, get_sample_name_and_extenstion, run

from reads_to_features_concat import reads_to_features_concat



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
        reads_to_features_concat(
            feature='exon',
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
) -> None:
    with open(new_tree_input, 'w') as f_out:
        for path in [*exon_concat_paths, starting_tree_input]:
            with file(path) as f_in:
                for line in f_in:
                    f_out.write(line)


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
    newick = fasta_to_newick(
        fasta_path=tree_input,
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
