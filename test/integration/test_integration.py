from importlib.metadata import metadata
import unittest

from src.transform import newick_to_imgs, reads_to_exons_concat, exons_concat_to_newick, reads_to_fastqc, alignment_to_flagstat
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

    def assertExpectedDirectoryFilesMatch(self, directory):
        exp_dir = join(EXP_DIR, directory)
        obs_dir = join(OBS_DIR, directory)
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

@unittest.skipIf(should_skip_integration_tests(), 'skipping integration test')
class IntegrationTestCase(unittest.TestCase, Assertions):
    
    def setUp(self):
        setUp()


class TestReadsToExonsConcat(IntegrationTestCase):

    def test_reads_to_exons_concat(self):
        reads_to_exons_concat(
            fastq=join(IN_DIR, 'isolate_1.fastq.gz'),
            reference='data/reference/pst-130_388_genes.fasta',
            gff='data/reference/pst-130_388_genes_as_positive_strand_landmarks.gff3',
            out_dir=join(OBS_DIR, 'isolate_1'),
        )
        self.assertExpectedDirectoryFilesMatch('isolate_1')

class TestExonsConcatToNewick(IntegrationTestCase):

    def test_exons_concat_to_newick(self):
        # TODO: this is a fragile way to compare trees, something with
        #       BIO.Phylo which re-roots, ladderizes, and does approximate
        #       comparison of branch lengths would safer.
        input_name = '6_isolates_8000_bases.fasta'
        exp_path = join(EXP_DIR, f'RAxML_bestTree.{splitext(input_name)[0]}.newick')
        obs_path = exons_concat_to_newick(join(IN_DIR, input_name), OBS_DIR)
        self.assertExpectedFilesMatch(exp_path, obs_path)

    def test_exons_concat_to_newick_exceptions(self):
        # When run with a fasta file, RAxML reports it cannot be parsed as a phylip.
        # This message can make it harder to find the actual error if RAxML crashes
        # so exons_concat_to_newick is supposed to supress that message.
        exception = ''
        try:
            exons_concat_to_newick(join(IN_DIR, 'mixed_length_sequences.fasta'), OBS_DIR)
        except Exception as e:
            exception = str(e)
        self.assertNotIn(
            'parse the alignment file as phylip file',
            exception,
            'warning about RAxML input file being fasta instead of phylip should have been suprressed but was not'
        )

class TestNewickToImgs(IntegrationTestCase):

    def test_newick_to_imgs(self):
        # TODO: this is a fragile way to compare plots, something with
        #       image histograms might be safer but its a tricky one.
        name = '50_isolates'
        newick_path = join(IN_DIR, f'RAxML_bestTree.{name}.newick')
        metadata_path = 'data/metadata_264_isolates.xlsx'

        newick_to_imgs(newick_path, metadata_path, join(OBS_DIR, f'{name}_imgs'), 'png')

        self.assertExpectedDirectoryFilesMatch(f'{name}_imgs')

class TestReadsToFastqc(IntegrationTestCase):

    def test_reads_to_fastqc(self):
        # TODO: could also test the zip file contents match, so that
        #       multiqc can pick up the data more easily, but it may
        #       be better to just test multiqc output instead
        obs_dir = join(OBS_DIR, 'fastqc')
        reads_to_fastqc(join(IN_DIR, 'isolate_1.fastq.gz'), obs_dir)
        obs_path = join(obs_dir, 'isolate_1_fastqc.html')
        exp_path = join(EXP_DIR, 'fastqc', 'isolate_1_fastqc.html')
        self.assertExpectedFilesMatch(exp_path, obs_path)

class TestAlignmentToFlagstat(IntegrationTestCase):

    def test_alignment_to_flagstat(self):
        alignment_to_flagstat(join(IN_DIR, 'isolate_1_sorted.bam'), join(OBS_DIR, 'flagstat'))
        self.assertExpectedDirectoryFilesMatch('flagstat')

if __name__ == '__main__':
    unittest.main()
