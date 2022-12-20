from os import makedirs
from os.path import isfile, join
from utils import get_sample_name_and_extenstion, run, file
from Bio.SeqIO import parse

def filter_fastq_read_length(
    fastq_in: str,
    fastq_out: str,
    max_read_length: int,
) -> None:
    with file(fastq_out, 'wt') as f_out:
        for record in parse(fastq_in, 'fastq'):
            if len(record.seq) <= max_read_length:
                f_out.write(record.format('fastq'))

# Make a pileup and return the path to it
def reads_to_pileup(
    fastq: str,
    reference: str,
    out_dir: str,
    threads=1,
    trim=True,
    max_read_length=None,
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

    if max_read_length:
        filtered = join(out_dir, f'{sample_name}_len_le_{max_read_length}{sample_ext}')
        print('filtering', end=' ', flush=True)
        filter_fastq_read_length(trimmed, filtered, max_read_length)
    else:
        filtered = trimmed

    print('aligning', end=' ', flush=True)
    # Don't recreate the index if it already exists as it may crash in parallel
    if not isfile(f'{reference}.bwt'):
        run(['bwa', 'index', reference])

    aligned_sam = join(out_dir, f'{sample_name}.sam')
    run(['bwa', 'mem', '-t', threads, reference, filtered], aligned_sam)
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
