#!/bin/bash

# Bash strict mode
set -euo pipefail

marple_pst_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

[ -d "$marple_pst_dir"/example ] && rm -r "$marple_pst_dir"/example
tar -xzf "$marple_pst_dir"/.example_backup.tar.gz 
