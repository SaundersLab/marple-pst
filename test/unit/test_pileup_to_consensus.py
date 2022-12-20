import unittest
from src.pileup_to_consensus import (
    base_ratios_from_reads,
    genotype_and_valid_base_ratios,
    pileup_row_to_consensus_at_locus
)

class TestPileupToConsensus(unittest.TestCase):

    def test_base_ratios_from_reads(self):
        cases = [
            ('G', 22, '.' * 22,           {'G': 1}),
            ('A', 10, '.' * 9 + 'T',      {'A': .9, 'T': .1}),
            ('A', 16, '....TTTTGGGGGGGG', {'A': .25, 'T': .25, 'G': .5}),
            ('T', 16, '....AAAAGGGGGGGG', {'A': .25, 'T': .25, 'G': .5}),
            ('A', 10, '.,.,.,TTTT',       {'A': .6, 'T': .4}),
            ('A', 10, '.,.,.,TTTT$',      {'A': .6, 'T': .4}),
            ('A', 10, '.,.,.,TTTT^',      {'A': .6, 'T': .4}),
            # Test indels removed
            ('A', 10, '.,.,.,+3ATATTTT',  {'A': .6, 'T': .4}),
            ('A', 10, '.,.,.,-2TCTTTT',   {'A': .6, 'T': .4}),
            # ^G. marks the start of a read segment where the ASCII 
            # minus 33 of G is the mapping quality and the character
            # AFTER those 2 is the actual read - in this case a refrence match
            ('T', 10, '.....^G.....',     {'T': 1.0}), 
        ]
        for reads, depth, ref, base_ratios in cases:
            self.assertEqual(
                base_ratios_from_reads(reads, depth, ref),
                base_ratios
            )

    def test_genotype_and_valid_base_ratios(self):
        cases = [
            # 0/0 : pref ≥ hetero_max or only refbase qualified
            ('A', {'A': 1}, '0/0', {'A': 1.0}),
            ('A', {'A': .9, 'T': .1}, '0/0', {'A': .9}),
            ('A', {'T': .1, 'A': .9}, '0/0', {'A': .9}),
            # 0/1 : ref base and 1 alt base qualified
            ('T', {'C': .217, 'T': .783}, '0/1', {'C': .217, 'T': .783}),
            ('A', {'A': .5, 'T': .5}, '0/1', {'A': .5, 'T': .5}),
            ('A', {'T': .5, 'A': .5}, '0/1', {'T': .5, 'A': .5}),
            # 1/1 : palt 1 ≥ hetero_max or only 1 alt base qualified
            ('T', {'T': .178, 'C': .822}, '1/1', {'C': .822}),
            ('A', {'A': .1, 'T': .9}, '1/1', {'T': .9}),
            ('T', {'A': .9, 'T': .1}, '1/1', {'A': .9}),
            ('A', {'A': .15, 'T': .15, 'G': .7}, '1/1', {'G': 0.7}),
            # 1/2 : 2 alt bases qualified and ref base did not
            ('C', {'A': .476, 'T': .524}, '1/2', {'A': .476, 'T': .524}),
            ('A', {'T': .5, 'C': .5}, '1/2', {'T': .5, 'C': .5}),
            # ? : 3+ bases qualified
            ('T', {'A': .2, 'C': .6, 'T': .2}, '?', {'A': .2, 'C': .6, 'T': .2}),
            ('A', {'A': .3, 'C': .4, 'T': .3}, '?', {'A': .3, 'C': .4, 'T': .3}),
            ('G', {b: .25 for b in 'ACTG'}, '?', {b: .25 for b in 'ACTG'}),
            # None : only ref has min <= freq <= max
            ('A', {'A': 0.7, 'T': 0.15, 'G': 0.15}, None, None),
        ]
        for ref, base_ratios, genotype, valid_base_ratios in cases:
            self.assertEqual(
                genotype_and_valid_base_ratios(ref, base_ratios, .2, .8),
                (genotype, valid_base_ratios)
            )

        # No bases qualified
        self.assertEqual(
            genotype_and_valid_base_ratios('G', {b: .25 for b in 'ACTG'}, .3, .7),
            (None, None)
        )

    def test_pileup_row_to_consensus(self):
        pileup_row_to_consensus_args = {
            'min_ref_depth': 2,
            'min_snp_depth': 10,
            'hetero_min': .2,
            'hetero_max': .8,
        }
        cases = [
            ('Chr1	1	A	10	.,.,.,TTTT	85.9..//.3', 'Chr1', 0, 'W'),
        ]
        for row, *consensus_at_locus in cases:
            self.assertEqual(
            pileup_row_to_consensus_at_locus(row, **pileup_row_to_consensus_args),
            tuple(consensus_at_locus)
        )
