import unittest

from src.transform import _base_ratios_from_reads, _get_genotype_and_valid_bases_and_valid_ratios

class TestTransform(unittest.TestCase):

    def test_base_ratios_from_reads(self):
        cases = [
            ('.' * 22,           22, 'G', {'G': 1}),
            ('.' * 9 + 'T',      10, 'A', {'A': .9, 'T': .1}),
            ('....TTTTGGGGGGGG', 16, 'A', {'A': .25, 'T': .25, 'G': .5}),
            ('....AAAAGGGGGGGG', 16, 'T', {'A': .25, 'T': .25, 'G': .5}),
            ('.,.,.,TTTT',       10, 'A', {'A': .6, 'T': .4}),
            ('.,.,.,TTTT$',      10, 'A', {'A': .6, 'T': .4}),
            ('.,.,.,TTTT^',      10, 'A', {'A': .6, 'T': .4}),
            ('.,.,.,+3ATATTTT',  10, 'A', {'A': .6, 'T': .4}),
            ('.,.,.,-2TCTTTT',   10, 'A', {'A': .6, 'T': .4}),
            # ^G. marks the start of a read segment where the ASCII 
            # minus 33 of G is the mapping quality and the character
            # AFTER those 2 is the actual read - in this case a refrence match
            ('.....^G.....',     10, 'T', {'T': 1.0}), 
        ]
        for reads, depth, ref, base_ratios in cases:
            self.assertEqual(_base_ratios_from_reads(reads, depth, ref), base_ratios)

    def test_get_genotype_and_valid_bases_and_valid_ratios(self):
            # 1/1	        palt 1 ≥ hetero_max or only 1 alt base qualified	T CC 1/1 0.822
            # 1/2	        2 alt bases qualified and ref base did not	        C AT 1/2 0.476,0.524
            # ?	        3+ bases qualified	                                T ACT ? 0.2,0.6,0.2
        cases = [
            # 0/0 : pref ≥ hetero_max or only refbase qualified
            ('A', ['A'],         ['1'],                 .2, .8, '0/0', 'AA',   '1'),
            ('A', ['A', 'T'],    ['0.9', '0.1'],        .2, .8, '0/0', 'AA',   '0.9'),
            ('A', ['T', 'A'],    ['0.1', '0.9'],        .2, .8, '0/0', 'AA',   '0.9'),
            # 0/1 : ref base and 1 alt base qualified
            ('T', ['C', 'T'],    ['0.217', '0.783'],    .2, .8, '0/1', 'CT',   '0.217,0.783'),
            ('A', ['A', 'T'],    ['0.5', '0.5'],        .2, .8, '0/1', 'AT',   '0.5,0.5'),
            ('A', ['T', 'A'],    ['0.5', '0.5'],        .2, .8, '0/1', 'TA',   '0.5,0.5'),
            # 1/1 : palt 1 ≥ hetero_max or only 1 alt base qualified
            ('T', ['T', 'C'],    ['0.178', '0.822'],    .2, .8, '1/1', 'CC',   '0.822'),
            ('A', ['A', 'T'],    ['0.1', '0.9'],        .2, .8, '1/1', 'TT',   '0.9'),
            ('T', ['A', 'T'],    ['0.9', '0.1'],        .2, .8, '1/1', 'AA',   '0.9'),
            # 1/2 : 2 alt bases qualified and ref base did not
            ('C', ['A', 'T'],    ['0.476', '0.524'],    .2, .8, '1/2', 'AT',   '0.476,0.524'),
            ('A', ['T', 'C'],    ['0.5', '0.5'],        .2, .8, '1/2', 'TC',   '0.5,0.5'),
            # ? : 3+ bases qualified
            ('T', ['A','C','T'], ['0.2', '0.6', '0.2'], .2, .8, '?',   'ACT',  '0.2,0.6,0.2'),
            ('A', ['A','C','T'], ['0.3', '0.4', '0.3'], .2, .8, '?',   'ACT',  '0.3,0.4,0.3'),
            ('G', list('ACTG'),  ['0.25'] * 4,          .2, .8, '?',   'ACTG', '0.25,0.25,0.25,0.25'),
        ]
        for ref, bases, ratio_strs, hetero_min, hetero_max, genotype, valid_bases, valid_ratios in cases:
            self.assertEqual(
                _get_genotype_and_valid_bases_and_valid_ratios(ref, bases, ratio_strs, hetero_min, hetero_max),
                (genotype, valid_bases, valid_ratios)
            )


if __name__ == '__main__':
    unittest.main()
