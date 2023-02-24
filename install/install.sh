#!/bin/bash

# Bash strict mode
set -euo pipefail

install_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
marple_pst_dir=$(dirname $install_dir)

if [ -d "$marple_pst_dir"/conda ] ; then
    while true; do
        read -p "WARNING: You may have already installed dependencies, do you want to overwrite them? (y/n)" yn
        case $yn in
            [Yy]* ) echo "INFO: Uninstalling and reinstalling" ; "$install_dir"/uninstall.sh ; break;;
            [Nn]* ) echo "INFO: Will not overwrite installation. Exiting." ; exit;;
            * ) ;;
        esac
    done
fi

echo "INFO: Installing into $marple_pst_dir"

full_operating_system_name="$(uname -s)"
case "${full_operating_system_name}" in
    Linux*)     operating_system_name=Linux;;
    Darwin*)    operating_system_name=Mac;;
    *)          operating_system_name="unsupported"
esac

if [ "$operating_system_name" = "unsupported" ]; then
    echo "ERROR: only Linux x86_64, macOS x86_64, and macOS arm64 are supported" >&2
    exit 1
fi

# try installing with curl, otherwise try installing with wget
which curl && curl -L "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" -o $install_dir/Mambaforge.sh
which curl || wget -O $install_dir/Mambaforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"

bash $install_dir/Mambaforge.sh -b -p "$marple_pst_dir/conda"
source "$marple_pst_dir/conda/etc/profile.d/conda.sh"
conda activate


if conda env list | grep '/marple-pst$' >/dev/null 2>&1 ; then 
    while true; do
        read -p "WARNING: You already have a conda environment called marple-pst, do you want to overwrite it? (y/n)" yn
        case $yn in
            [Yy]* ) echo "INFO: Overwriting existing marple-pst envirionment" ; break;;
            [Nn]* ) echo "INFO: Will not overwrite existing marple-pst envirionment. Exiting." ; exit;;
            * ) ;;
        esac
    done
fi

machine_name="$(uname -m)"
if [ "$operating_system_name" = "Mac" ] && [ "$machine_name" = "arm64" ]; then
    CONDA_SUBDIR=osx-64 mamba env create --force -f env.yml
    conda activate marple-pst
    conda config --env --set subdir osx-64
else 
    mamba env create --force -f env.yml
    conda activate marple-pst
fi

./test/run_tests.sh

