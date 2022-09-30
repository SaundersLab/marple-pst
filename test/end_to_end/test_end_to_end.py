import unittest
from subprocess import run
import filecmp
import os
from os.path import basename, join

def should_skip_end_to_end_tests():
    try:
        return os.environ['SKIP_END_TO_END'] == 'true'
    except:
        return False

def matches(s, matches=[], suffixes=[], prefixes=[]):
    if s in matches:
        return True
    if any([s.endswith(suffix) for suffix in suffixes]):
        return True
    if any([s.startswith(prefix) for prefix in prefixes]):
        return True
    return False

def should_ignore_file(file):
    return matches(
        file,
        matches=['.DS_Store'],
        suffixes=['.bam', '.pdf', '.html', 'xlsx', '.zip'],
        prefixes=['RAxML_parsimonyTree', 'RAxML_info', 'RAxML_log', 'RAxML_result'],
    )

def should_ignore_directory(directory):
    return matches(
        directory,
        matches=['multiqc_data', 'fastq'],
        prefixes=['barcode'],
    )

class Assertions:

    def assertFileExists(self, obs_path):
        if not os.path.isfile(obs_path):
            raise AssertionError(f'expected output file {obs_path} was not created')

    def assertExpectedFilesMatch(self, exp_path, obs_path):
        exp_name = basename(exp_path)
        self.assertFileExists(obs_path)
        if not filecmp.cmp(exp_path, obs_path):
            raise AssertionError(f'contents of {exp_name} does not match the expected output')

@unittest.skipIf(should_skip_end_to_end_tests(), 'skipping end to end test')
class TestRunExample(unittest.TestCase, Assertions):

    def test_run_example(self):
        run('./run_example.sh')
        for exp_dir, _, files in os.walk('finished_example'):
            if should_ignore_directory(basename(exp_dir)):
                continue
            obs_dir = join('example', exp_dir[len('finished_example') + 1:])
            files = [file for file in files if not should_ignore_file(file)]
            for file in files:
                self.assertExpectedFilesMatch(
                    join(exp_dir, file),
                    join(obs_dir, file),
                )

if __name__ == '__main__':
    unittest.main()
