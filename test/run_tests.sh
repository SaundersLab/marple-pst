#!/bin/bash

flag=$1

[ "$flag" == "SKIP_INTEGRATION" ] && export SKIP_INTEGRATION=true 

[ "$SKIP_INTEGRATION" == true ] && echo "skiping integtation tests"

export PYTHONPATH=$PYTHONPATH:$PWD/src

[ "$flag" != "COVERAGE" ] && python3 -m unittest

[ "$flag" == "COVERAGE" ] && coverage run --omit=src/__init__.py --source=src -m unittest discover && coverage report -m

exit 0
