from utils import get_sample_name_and_extenstion, file, write_fasta
from collections import defaultdict
from os.path import join
from Bio.SeqIO import parse


def consensus_to_exons(
    consensus_path: str,
    gff: str,
    exons_path: str,
) -> None:
    """
     Extract the exons from the gene sequence. Take shortcuts because
     we have a very specific GFF file. The GFF and reference were
     transformed so everything is on plus strand, and seqid is the gene.
     So we don't have to worry about phase or attributes.
    """
    # Find the exon positions
    exon_positions = defaultdict(set)
    with file(gff) as f:
        for line in f:
            # Assumes everything's on the positive strand
            seqid, _, type_, start, end, _, strand, _, _ = line.strip().split('\t')
            assert strand == '+', 'Only + strand supported'
            if type_ == 'exon':
                exon_positions[seqid].update(
                    set(range(int(start) - 1, int(end))))

    # Write the gene bases only in the exon positions
    consensus_exons = {
        r.id: ''.join(r.seq[i] for i in sorted(exon_positions[r.id]))
        for r in parse(consensus_path, 'fasta')
    }
    
    write_fasta(consensus_exons, exons_path, sort=True)
