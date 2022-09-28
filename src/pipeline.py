from transform import reads_to_exons_concat, sample_report
from asyncio import streams
import sys
from os.path import join, dirname, realpath
from typing import List
import os
from utils import run

def pipeline(
    fastq_paths: List[str],
    reference: str,
    gff: str,
    out_dirs: List[str],
    multiqc_config: str,
):
    for fastq, out_dir in zip(fastq_paths, out_dirs):
        print(fastq)
        reads_to_exons_concat(
            fastq=fastq,
            reference=reference,
            gff=gff,
            out_dir=out_dir,
        )
        sample_report(out_dir)
    report_dirs = [join(out_dir, 'report') for out_dir in out_dirs]
    run(['multiqc', '--config', multiqc_config] + report_dirs, out='/dev/null')


if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    reference_dir = join(marple_dir, 'data', 'reference')
    reference = join(reference_dir, 'pst-130_388_genes.fasta')
    gff = join(reference_dir, 'pst-130_388_genes_as_positive_strand_landmarks.gff3')
    relative_fastq_paths = sys.argv[1:]
    fastq_paths = [realpath(path) for path in relative_fastq_paths]
    out_dirs = [dirname(path) for path in fastq_paths]
    multiqc_config = join(marple_dir, 'config', 'multiqc_config.yaml')
    pipeline(fastq_paths, reference, gff, out_dirs, multiqc_config)

