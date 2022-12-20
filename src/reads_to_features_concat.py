from reads_to_pileup import reads_to_pileup
from utils import get_sample_name_and_extenstion
from os.path import join
from pileup_to_consensus import pileup_to_consensus
from consensus_to_features import consensus_to_features

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