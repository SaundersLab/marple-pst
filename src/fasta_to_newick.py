import tempfile
from os import makedirs
from os.path import abspath, join
from utils import get_sample_name_and_extenstion, pushd, run

def fasta_to_newick(exons_concat: str, out_dir: str, n_threads=1) -> str:
    makedirs(out_dir, exist_ok=True)
    collection_name, _ = get_sample_name_and_extenstion(exons_concat, 'fasta')
    # need absolute path because we're about to run raxml from the output directory
    exons_concat = abspath(exons_concat)

    # output name cannot contain slashes, so run raxml from the output directory
    with pushd(out_dir):
        # make a temporary file to prevent raxml log being written to screen
        with tempfile.TemporaryFile() as log:
            try:
                args = [
                    'raxmlHPC-PTHREADS-SSE3',
                     '-s', exons_concat,
                     '-m', 'GTRGAMMA',
                     '-n', f'{collection_name}.newick',
                     '-p', '100',
                ]
                if n_threads >= 2:
                    args += ['-T', str(n_threads)]
                run(args, out=log)
            except:
                # if raxml failed then raise an exception with the error message
                # which was written to the temporary file (redirected from stdout)
                log.seek(0)
                error = log.read().decode()
                # When run with a fasta file, RAxML reports it cannot be parsed as a phylip.
                # This message can make it harder to find the actual error if RAxML crashes,
                # so supress that message.
                error = error.replace(
                    "\nRAxML can't, parse the alignment file as phylip file \nit will now try to parse it as FASTA file\n\n",
                    ''
                )
                raise Exception(error)

    return join(out_dir, f'RAxML_bestTree.{collection_name}.newick')

