import colorsys
import contextlib
import gzip
import subprocess
from hashlib import sha1
from os import chdir, getcwd
from os.path import basename
from typing import Dict, IO, List, Tuple, Union

import matplotlib


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


def run(args: List[str], out: Union[str, IO] = None) -> None:
    if out is not None:
        if isinstance(out, str):
            with open(out, 'wt') as f:
                process = subprocess.run(
                    args, text=True, stdout=f, stderr=subprocess.PIPE)
        else:
            process = subprocess.run(
                args, text=True, stdout=out, stderr=subprocess.PIPE)
    else:
        process = subprocess.run(args, text=True, stderr=subprocess.PIPE)
    if process.returncode:
        raise Exception(process.stderr)


def get_file_extenstion(path: str, candidate_exts: List[str]) -> str:
    for ext in candidate_exts:
        if path.endswith(ext):
            return ext
    extensions_str = '"' + '", "'.join(candidate_exts) + '"'
    raise ValueError(f'Unknown extension for file "{path}". Tried these extenstions: {extensions_str}')


def get_sample_name_and_extenstion(path: str, candidate_exts: Union[str, List[str]]) -> Tuple[str, str]:
    if isinstance(candidate_exts, str):
        candidate_exts = {
            'fastq': ['.fastq.gz', '.fq.gz', '.fastq', '.fq'],
            'pileup': ['.pileup.gz', '.mpileup.gz', '.pileup', '.mpileup'],
            'fasta': ['.fa.gz', '.fasta.gz', '.fna.gz', '.ffn.gz', '.fa', '.fasta', '.fna', '.ffn'],
            'alignment': ['.bam', '.sam', '.cram'],
            'newick': ['.newick.gz', '.tree.gz', '.newick', '.tree'],
            'snp_ratios': ['_snp_ratios.tsv'],
        }[candidate_exts]
    sample_filename = basename(path)
    sample_ext = get_file_extenstion(sample_filename, candidate_exts)
    sample_name = sample_filename[:-len(sample_ext)]
    raxml_prefix = 'RAxML_bestTree.'
    if sample_name.startswith(raxml_prefix):
        sample_name = sample_name[len(raxml_prefix):]
    return sample_name, sample_ext


@contextlib.contextmanager
def pushd(new_dir: str):
    previous_dir = getcwd()
    chdir(new_dir)
    try:
        yield
    finally:
        chdir(previous_dir)


def string_to_color(s: str) -> str:
    if s == '?':
        return '#00DD00'
    hash_ = str(int(sha1(str(s).encode("utf-8")).hexdigest(), 16))
    color_tuple = tuple(
        (int(hash_[(i * 3):(i * 3) + 3]) % 255) / 300 for i in range(3))
    return matplotlib.colors.to_hex(color_tuple)


def darken_color(color_hex: str) -> str:
    h, l, s = matplotlib.colors.hex2color(color_hex)
    return matplotlib.colors.to_hex(colorsys.hls_to_rgb(h, .25, s))

def write_fasta(fasta: Dict[str, str], path: str) -> None:
    with file(path, 'wt') as f_out:
        for header in sorted(fasta):
            assert isinstance(header, str)
            assert isinstance(fasta[header], str)
            f_out.write('>' + header + '\n' + fasta[header] + '\n')
