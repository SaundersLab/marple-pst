import unittest

from src.transform import reads_to_exons_concat
from shutil import rmtree
import filecmp
import os

def should_skip_integration_tests():
    try:
        return os.environ['SKIP_INTEGRATION'] == 'true'
    except:
        return False


class Assertions:

    def assertExpectedDirectoryFilesMatch(self, exp_dir, obs_dir):
        dircmp = filecmp.dircmp(exp_dir, obs_dir)
        missing = dircmp.left_only
        if missing:
            raise AssertionError(f'{len(missing)} expected output file(s) were not created: ' + "\n".join(missing))
        match, mismatch, errors = filecmp.cmpfiles(exp_dir, obs_dir, dircmp.left_list, False)
        if mismatch:
            raise AssertionError(f'{len(mismatch)} file(s) do not match the expected output: ' + "\n".join(mismatch))
        if errors:
            raise AssertionError(f'{len(errors)} file(s) could not be compared: ' + "\n".join(errors))


class TestReadsToExonsConcat(unittest.TestCase, Assertions):

    def setUp(self):
        rmtree('test/integration/observed', ignore_errors=True)

    @unittest.skipIf(should_skip_integration_tests(), 'skipping integration test')
    def test_reads_to_exons_concat(self):

        exp_dir = 'test/integration/expected/isolate_1'
        obs_dir = 'test/integration/observed/isolate_1'

        reads_to_exons_concat(
            fastq='test/integration/input/isolate_1.fastq.gz',
            reference='data/reference/pst-130_388_genes.fasta',
            gff='data/reference/pst-130_388_genes_as_positive_strand_landmarks.gff3',
            out_dir=obs_dir,
        )

        self.assertExpectedDirectoryFilesMatch(exp_dir, obs_dir)


if __name__ == '__main__':
    unittest.main()
