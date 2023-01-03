from os.path import dirname, join, realpath
import argparse
from typing import List, Dict
from utils import file
from fasta_to_newick import fasta_to_newick
from os import makedirs
from newick_to_images import newick_to_images

def exon_concat_paths_to_tree_input(
    exon_concat_paths: List[str],
    starting_tree_input: str,
    new_tree_input: str,
) -> None:
    with open(new_tree_input, 'w') as f_out:
        for path in [*exon_concat_paths, starting_tree_input]:
            with file(path) as f_in:
                for line in f_in:
                    f_out.write(line)


def exon_concat_paths_to_tree_imgs(
    exon_concat_paths: List[str],
    starting_tree_input: str,
    tree_name: str,
    out_dir: str,
    metadata_path: str,
    n_threads=1,
    img_fmt='pdf',
) -> Dict[str, str]:
    makedirs(out_dir, exist_ok=True)
    tree_input = join(out_dir, f'{tree_name}.fasta')
    exon_concat_paths_to_tree_input(
        exon_concat_paths=exon_concat_paths,
        starting_tree_input=starting_tree_input,
        new_tree_input=tree_input,
    )
    print('Tree: making', end=' ', flush=True)
    newick = fasta_to_newick(
        fasta_path=tree_input,
        out_dir=out_dir,
        n_threads=n_threads,
    )
    print('visualising', flush=True)
    imgs = newick_to_images(
        newick_path=newick,
        metadata_path=metadata_path,
        out_dir=out_dir,
        img_fmt=img_fmt,
    )
    return imgs


if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')

    parser = argparse.ArgumentParser(description='Create and visualise tree with new samples')
    parser.add_argument(
        'relative_exon_concat_paths',
        type=str,
        nargs='*',
        help='Concatenated exons file for new samples to add to the tree',
        default=[],
    )
    parser.add_argument(
        '--meta',
        help='Path to spreadsheet containing isolate metadata and styles',
        default=join(data_dir, 'metadata_264_isolates.xlsx')
    )
    parser.add_argument(
        '--start',
        help='Starting tree input, i.e. concatenated exons file for previously sequenced samples',
        default=join(data_dir, '57_isolates_388_genes_exons.fasta.gz')
    )
    parser.add_argument('--name', help='Name to use as a prefix for all created tree files',)
    parser.add_argument('--out_dir', help='Directory to create tree files in', default=realpath('.'))
    parser.add_argument('--threads', type=int, help='Number of threads to use for RAxML', default=1)
    parser.add_argument('--img_fmt', help='Format to output tree images as', default='pdf')

    args = parser.parse_args()
    exon_concat_paths_to_tree_imgs(
        exon_concat_paths=[realpath(path) for path in args.relative_exon_concat_paths],
        starting_tree_input=args.start,
        tree_name=args.name,
        out_dir=args.out_dir,
        metadata_path=args.meta,
        n_threads=args.threads,
        img_fmt=args.img_fmt,
    )
