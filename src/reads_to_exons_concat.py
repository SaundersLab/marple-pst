from os.path import join, dirname, realpath
import argparse
from utils import get_sample_name_and_extenstion, run
from typing import List
from reads_to_features_concat import reads_to_features_concat
from report import report

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


if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')
    reference_dir = join(data_dir, 'reference')

    parser = argparse.ArgumentParser(description='Align reads to reference genes then extract and concatenate exons')
    parser.add_argument(
        'relative_fastq_paths',
        type=str,
        nargs='+',
        help='Read files for aligning to the reference genes. Output for each read will be created in the same directory as the fastq file',
    )
    parser.add_argument(
        '--ref',
        help='Reference FASTA file to align reads to',
        default=join(reference_dir, 'pst-130_388_genes.fasta'),
    )
    parser.add_argument(
        '--gff',
        help='Annotation file giving the position of exons within the reference genes. Genes should be specificied as if everything was on the + strand.',
        default=join(reference_dir, 'pst-130_388_genes_as_positive_strand_landmarks.gff3'),
    )
    parser.add_argument(
        '--multiqc_config',
        help='MultiQC config file for creating the report',
        default=join(marple_dir, 'config', 'multiqc_config.yaml'),
    )
    parser.add_argument('--threads', type=int, help='Number of threads to use', default=2)
    parser.add_argument('--trim', help='Should FASTQ files be trimmed (yes/no)', default='yes')
    parser.add_argument('--max_read_length', help='Maximum length of read to align', default=None)
    args = parser.parse_args()
    fastq_paths = [realpath(path) for path in args.relative_fastq_paths]
    out_dirs = [dirname(path) for path in fastq_paths]
    reads_list_to_exons_concat_with_report(
        fastq_paths=fastq_paths,
        reference=args.ref,
        gff=args.gff,
        out_dirs=out_dirs,
        multiqc_config=args.multiqc_config,
        threads=args.threads,
        trim=args.trim.lower().startswith('y'),
    )
