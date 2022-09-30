from os.path import dirname, join, realpath
from transform import exon_concat_paths_to_tree_imgs
import argparse



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
        default=join(data_dir, '56_isolates_388_genes_exons.fasta.gz')
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
