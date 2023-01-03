from os import makedirs
from os.path import join, isfile
from typing import Tuple, Iterable
from utils import run, get_sample_name_and_extenstion
from Bio.SeqIO import parse
import pandas as pd

def reads_to_fastqc(fastq: str, out_dir: str) -> Tuple[str, str]:
    makedirs(out_dir, exist_ok=True)
    sample_name, _ = get_sample_name_and_extenstion(fastq, 'fastq')
    fastqc_page = join(out_dir, f'{sample_name}_fastqc.html')
    fastqc_data = join(out_dir, f'{sample_name}_fastqc.zip')
    run(['fastqc', '--quiet', '-o', out_dir, fastq])
    return fastqc_page, fastqc_data


def alignment_to_flagstat(alignment: str, flagstat_path: str) -> None:
    run(['samtools', 'flagstat', alignment], flagstat_path)

def seq_to_pct_coverage(seq: str) -> float:
    seq = str(seq).upper()
    n_unknown = seq.count('?') + seq.count('N')
    pct_unknown = 100 * n_unknown / len(seq)
    return 100 - pct_unknown

def consensus_to_coverage(consensus_path: str, coverage_path: str, step=1) -> None:
    cov_pcts = [seq_to_pct_coverage(r.seq) for r in parse(consensus_path, 'fasta')]

    pct_coverage_thresholds = list(range(0, 101, step))
    number_of_genes_with_pct_coverage_ge_thresholds = [
        len([pct for pct in cov_pcts if pct >= min_pct])
        for min_pct in pct_coverage_thresholds
    ]
    pd.DataFrame({
        'pct_exon_positions_covered': pct_coverage_thresholds,
        'number_of_genes': number_of_genes_with_pct_coverage_ge_thresholds,
    }).to_csv(coverage_path, index=None, header=None)

# This is a work around to ensure that multiqc can create the plot specified in
# config/multiqc_config.yaml.
# Unless there is at least one _mqc.yaml in the directory, multiqc will not 
# pick up any of the custom_data specified in the --config file. To ensure it
# does pick up those files, we create a file for each sample with only id.
def create_empty_config_required_for_gene_coverage_mqc(sample: str, out_dir: str):
    config_name = f'{sample}_empty_config_required_for_gene_coverage_mqc'
    with open(join(out_dir, f'{config_name}.yaml'), 'w') as f:
        f.write(f'id: "{config_name}"' + '\n')


def report(sample_dir: str, sample_name: str):
    out_dir = join(sample_dir, 'report')
    for fastq_ext in ['.fastq', '.fastq.gz', '.fq', '.fq.gz']:
        fastq_path = join(sample_dir, f'{sample_name}{fastq_ext}')
        if isfile(fastq_path):
            reads_to_fastqc(fastq_path, out_dir)
            break
    flagstat_path = join(out_dir, f'{sample_name}.txt')
    alignment_to_flagstat(join(sample_dir, f'{sample_name}.bam'), flagstat_path)
    coverage_path = join(out_dir, f'{sample_name}.csv')
    consensus_to_coverage(join(sample_dir, f'{sample_name}.fasta'), coverage_path)
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
