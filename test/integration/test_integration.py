from importlib.metadata import metadata
import unittest

from src.transform import newick_to_imgs, reads_to_exons_concat, exons_concat_to_newick
from shutil import rmtree
import filecmp
import os
from os.path import basename, join, splitext

def should_skip_integration_tests():
    try:
        return os.environ['SKIP_INTEGRATION'] == 'true'
    except:
        return False

OBS_DIR = 'test/integration/observed'
EXP_DIR = 'test/integration/expected'
IN_DIR = 'test/integration/input'

def setUp():
    rmtree(OBS_DIR, ignore_errors=True)

class Assertions:

    def assertExpectedDirectoryFilesMatch(self, exp_dir, obs_dir):
        dircmp = filecmp.dircmp(exp_dir, obs_dir)
        missing = dircmp.left_only
        if missing:
            raise AssertionError(f'{len(missing)} expected output file(s) were not created: ' + "\n".join(missing))
        match, mismatch, errors = filecmp.cmpfiles(exp_dir, obs_dir, dircmp.left_list, False)
        if mismatch:
            raise AssertionError(f'contents of {len(mismatch)} file(s) do not match the expected outputs: ' + "\n".join(mismatch))
        if errors:
            raise AssertionError(f'{len(errors)} file(s) could not be compared: ' + "\n".join(errors))

    def assertExpectedFilesMatch(self, exp_path, obs_path):
        exp_name = basename(exp_path)
        if not os.path.isfile(obs_path):
            raise AssertionError(f'expected output file {exp_name} was not created')
        if not filecmp.cmp(exp_path, obs_path):
            raise AssertionError(f'contents of {exp_name} does not match the expected output')


class TestReadsToExonsConcat(unittest.TestCase, Assertions):

    def setUp(self):
        setUp()

    @unittest.skipIf(should_skip_integration_tests(), 'skipping integration test')
    def test_reads_to_exons_concat(self):

        exp_dir = join(EXP_DIR, 'isolate_1')
        obs_dir = join(OBS_DIR, 'isolate_1')

        reads_to_exons_concat(
            fastq=join(IN_DIR, 'isolate_1.fastq.gz'),
            reference='data/reference/pst-130_388_genes.fasta',
            gff='data/reference/pst-130_388_genes_as_positive_strand_landmarks.gff3',
            out_dir=obs_dir,
        )

        self.assertExpectedDirectoryFilesMatch(exp_dir, obs_dir)

class TestExonsConcatToNewick(unittest.TestCase, Assertions):

    def setUp(self):
        setUp()

    @unittest.skipIf(should_skip_integration_tests(), 'skipping integration test')
    def test_exons_concat_to_newick(self):
        # TODO: this is a fragile way to compare trees, something with
        #       BIO.Phylo which re-roots, ladderizes, and does approximate
        #       comparison of branch lengths would safer.
        input_name = '6_isolates_8000_bases.fasta'
        exp_path = join(EXP_DIR, f'RAxML_bestTree.{splitext(input_name)[0]}.newick')
        obs_path = exons_concat_to_newick(join(IN_DIR, input_name), OBS_DIR)
        self.assertExpectedFilesMatch(exp_path, obs_path)

class TestNewickToImgs(unittest.TestCase, Assertions):

    def setUp(self):
        setUp()

    @unittest.skipIf(should_skip_integration_tests(), 'skipping integration test')
    def test_newick_to_imgs(self):
        # TODO: this is a fragile way to compare plots, something with
        #       image histograms might be safer but its a tricky one.
        name = '50_isolates'
        newick_path = join(IN_DIR, f'RAxML_bestTree.{name}.newick')
        metadata_path = 'data/metadata_264_isolates.xlsx'

        obs_dir = join(OBS_DIR, f'{name}_imgs')
        newick_to_imgs(newick_path, metadata_path, obs_dir, 'png')

        exp_dir = join(EXP_DIR, f'{name}_imgs')
        self.assertExpectedDirectoryFilesMatch(exp_dir, obs_dir)



if __name__ == '__main__':
    unittest.main()
