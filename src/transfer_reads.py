

import sys
import glob
import os
import subprocess
import shutil

def du(path):
    return int(subprocess.check_output(['du','-s', path]).split()[0].decode('utf-8'))

marple_data_dir = "/mnt/c/Users/MARPLE/Documents/MARPLE_data"

experiment=sys.argv[1]
barcodes = {
    s.split('=')[0]: s.split('=')[1]
    for s in sys.argv[2:]
}

barcode_dirs = glob.glob(marple_data_dir + "/" + experiment + "/*/*/basecalling/pass/barcode*")
for barcode_dir in barcode_dirs:
    barcode = os.path.basename(barcode_dir)
    dir_size = du(barcode_dir)
    contains_reads = dir_size > 5000
    if not contains_reads:
        if barcode in barcodes:
            raise ValueError(f'An alias was given for {barcode} but could not find enough FASTQ data in {barcode_dir}')
        continue
    if barcode not in barcodes:
        raise ValueError(f'No alias given for {barcode} but it appears to have read data in {barcode_dir}')

for barcode_dir in barcode_dirs:
    barcode = os.path.basename(barcode_dir)
    if barcode not in barcodes:
        continue
    sample = barcodes[barcode]
    out_dir = experiment + '/' + sample
    os.makedirs(out_dir, exist_ok=True)
    out_path = out_dir + '/' + sample + '.fastq'
    with open(out_path, 'wb') as f_out:
        for fastq_path in sorted(glob.glob(barcode_dir + '/*.fastq')):
            with open(fastq_path, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)
