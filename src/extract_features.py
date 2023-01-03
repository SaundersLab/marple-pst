from utils import run, unwrap_fasta
from tempfile import NamedTemporaryFile

def extract_features(
    genomic_fasta_path: str,
    gff: str,
    features_path: str,
    feature: str,
) -> None:
    """
    Extract the features (exons or cds) from the gene sequence
    """
    try:
        gffread_flag = {'cds': '-x', 'exon': '-w'}[feature]
    except:
        raise ValueError(f'unsupported feature: {feature}')

    with NamedTemporaryFile() as wrapped_temp_file:
        run(['gffread', gffread_flag, wrapped_temp_file.name, '-g', genomic_fasta_path, gff])
        unwrap_fasta(wrapped_temp_file.name, features_path)
