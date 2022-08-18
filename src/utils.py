import gzip
from typing import IO, Tuple
import subprocess
from typing import List, Union
from os.path import basename

def file(path, mode='rt') -> IO:
    """Create a file object from path. Works regardless of compression based on extension.
    Parameters
    ----------
    path : str
        Path to file to open
    mode : str, default='rt'
        Any mode used in open or gzip.open
    Returns
    ----------
    opened_file : file object
    Examples
    --------
    >>> ','.join(file('f.txt.gz'))
    'a,b,c'
    >>> ','.join(file('f.txt'))
    'a,b,c'
    """
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def run(args: List[str], out: str = None) -> None:
    if out is not None:
        # subprocess stdout will not be compressed even if giving gzip.open
        with open(out, 'wt') as f:
            process = subprocess.run(args, text=True, stdout=f, stderr=subprocess.PIPE)
    else:
        process = subprocess.run(args, text=True, stderr=subprocess.PIPE)
    if process.returncode:
        raise Exception(process.stderr)

def get_file_extenstion(path: str, candidate_exts: List[str]) -> str:
    for ext in candidate_exts:
        if path.endswith(ext):
            return ext
    raise ValueError('Unknown extension for file ' + path)

def get_sample_name_and_extenstion(path: str, candidate_exts: Union[str, List[str]]) -> Tuple[str, str]:
    if isinstance(candidate_exts, str):
        candidate_exts = {
            'fastq': ['.fastq.gz', '.fq.gz', '.fastq', '.fq'],
            'pileup': ['.pileup.gz', '.mpileup.gz', '.pileup', '.mpileup'],
            'fasta': ['.fa.gz', '.fasta.gz', '.fna.gz', '.fa', '.fasta', '.fna'],
        }[candidate_exts]
    sample_filename = basename(path)
    sample_ext = get_file_extenstion(sample_filename, candidate_exts)
    sample_name = sample_filename[:-len(sample_ext)]
    return sample_name, sample_ext

