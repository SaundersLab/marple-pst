from utils import file, write_fasta
from collections import defaultdict
from Bio.SeqIO import parse


def extract_features(
    genomic_fasta_path: str,
    gff: str,
    features_path: str,
    feature: str,
) -> None:
    """
    Extract the features (exons or cds) from the gene sequence. Take shortcuts
    because we have a very specific GFF file. The GFF and reference were
    transformed so everything is on plus strand, and seqid is the gene,
    so we don't have to worry about phase or attributes.
    """
    assert feature in {'exon', 'cds'}, f'unsupported feature: {feature}'

    # Find the feature positions
    feature_positions = defaultdict(set)
    with file(gff) as f:
        for line in f:
            # Assumes everything's on the positive strand
            seqid, _, type_, start, end, _, strand, _, _ = line.strip().split('\t')
            assert strand == '+', 'Only + strand supported'
            if type_ == feature:
                feature_positions[seqid].update(
                    set(range(int(start) - 1, int(end)))
                )

    # Write the gene bases only in the feature positions
    features = {
        r.id: ''.join(r.seq[i] for i in sorted(feature_positions[r.id]))
        for r in parse(genomic_fasta_path, 'fasta')
    }

    write_fasta(features, features_path, sort=True)
