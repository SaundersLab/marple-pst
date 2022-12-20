import sys
from os.path import dirname, join, realpath
from newick_to_images import newick_to_images
import argparse

if __name__ == '__main__':
    src_dir = dirname(realpath(__file__))
    marple_dir = dirname(src_dir)
    data_dir = join(marple_dir, 'data')

    parser = argparse.ArgumentParser(description='Create and visualise tree')
    parser.add_argument(
        'newick_path',
        type=str,
        help='Tree in Newick format',
    )
    parser.add_argument(
        '--meta',
        help='Path to spreadsheet containing isolate metadata and styles',
        default=join(data_dir, 'metadata_264_isolates.xlsx')
    )
    parser.add_argument('--out_dir', help='Directory to create tree files in', default=None)
    parser.add_argument('--img_fmt', help='Format to output tree images as', default='pdf')
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = dirname(realpath(args.newick_path))
        print(args.out_dir)
    newick_to_images(
        newick_path=args.newick_path,
        metadata_path=args.meta,
        out_dir=args.out_dir,
        img_fmt=args.img_fmt,
    )
