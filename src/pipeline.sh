#!/bin/bash

src_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pst_dir=$( dirname "$src_dir" )
source "$marple_pst_dir"/marple_pst_miniconda/bin/activate marple-pst
conda activate marple-pst
python3 "$src_dir"/pipeline.py $@
