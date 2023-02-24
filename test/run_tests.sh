#!/bin/bash

export PYTHONPATH=$PYTHONPATH:$PWD/src

# Bash strict mode
set -euo pipefail

source conda/bin/activate marple-pst

coverage=false
# allow user to skip certain tests and get coverage if requested
while [[ $# -gt 0 ]] ; do
    flag=$1
    [ "$flag" == "skip_integration" ] && export SKIP_INTEGRATION=true 
    [ "$flag" == "skip_end_to_end" ] && export SKIP_END_TO_END=true 
    [ "$flag" == "coverage" ] && coverage=true
    shift 1
done


[ "$coverage" == "false" ] && python3 -m unittest
[ "$coverage" == "true" ] && coverage run --omit=src/__init__.py --source=src -m unittest discover && coverage report -m

exit 0
