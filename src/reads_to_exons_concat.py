from transform import reads_list_to_exons_concat_with_report
import sys
from os.path import join, dirname, realpath

if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    reference_dir = join(marple_dir, 'data', 'reference')
    reference = join(reference_dir, 'pst-130_388_genes.fasta')
    gff = join(reference_dir, 'pst-130_388_genes_as_positive_strand_landmarks.gff3')
    multiqc_config = join(marple_dir, 'config', 'multiqc_config.yaml')
    relative_fastq_paths = sys.argv[1:]
    fastq_paths = [realpath(path) for path in relative_fastq_paths]
    out_dirs = [dirname(path) for path in fastq_paths]
    reads_list_to_exons_concat_with_report(
        fastq_paths=fastq_paths,
        reference=reference,
        gff=gff,
        out_dirs=out_dirs,
        multiqc_config=multiqc_config,
    )
