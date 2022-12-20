from os import makedirs
from os.path import isfile, join
from utils import get_sample_name_and_extenstion, run

# Make a pileup and return the path to it
def reads_to_pileup(
    fastq: str,
    reference: str,
    out_dir: str,
    threads=1,
    trim=True,
) -> str:
    makedirs(out_dir, exist_ok=True)
    sample_name, sample_ext = get_sample_name_and_extenstion(fastq, 'fastq')
    threads = str(threads)
    if trim:
        trimmed = join(out_dir, f'{sample_name}_porechopped{sample_ext}')
        print('trimming', end=' ', flush=True)
        run(['porechop', '--threads', threads, '-i', fastq], trimmed)
    else:
        trimmed = fastq
    print('aligning', end=' ', flush=True)
    # Don't recreate the index if it already exists as it may crash in parallel
    if not isfile(f'{reference}.bwt'):
        run(['bwa', 'index', reference])
    aligned_sam = join(out_dir, f'{sample_name}.sam')
    run(['bwa', 'mem', '-t', threads, reference, trimmed], aligned_sam)
    aligned_bam = join(out_dir, f'{sample_name}_unsorted.bam')
    run(['samtools', 'view', '-@', threads, '-S', '-b', aligned_sam], aligned_bam)
    sorted_bam = join(out_dir, f'{sample_name}.bam')
    run(['samtools', 'sort', '-@', threads, aligned_bam], sorted_bam)
    # Don't recreate the index if it already exists as it may crash in parallel
    if not isfile(f'{reference}.fai'):
        run(['samtools', 'faidx', reference])
    pileup = join(out_dir, f'{sample_name}.pileup')
    run(['samtools', 'mpileup', '-f', reference, sorted_bam], pileup)

    # Let the caller know where to find the pileup file
    return pileup
