from reads_to_pileup import reads_to_pileup
from utils import get_sample_name_and_extenstion, write_fasta
from os.path import join
from pileup_to_consensus import pileup_to_consensus
from extract_features import extract_features
from Bio.SeqIO import parse

def concat_fasta_sequences(
    fasta_path: str, concat_fasta_path: str, header: str
) -> None:
    """
    Concatenate all the sequences in a fasta file into a single record.
    """
    sequences_combined = ''.join(str(r.seq) for r in parse(fasta_path, 'fasta'))
    write_fasta({header: sequences_combined}, concat_fasta_path)

def reads_to_features_concat(
    feature: str,
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
    assert feature in {'exon', 'cds'}, f'unsupported feature: {feature}'
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
    plural = 'exons' if feature == 'exon' else 'cds'
    features_path = join(out_dir, f'{sample_name}_{plural}.fasta')
    extract_features(consensus_path, gff, features_path, feature=feature)
    features_concat_path = join(out_dir, f'{sample_name}_{plural}_concat.fasta')
    concat_fasta_sequences(features_path, features_concat_path, sample_name)
    return features_concat_path
