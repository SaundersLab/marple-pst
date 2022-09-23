#!/bin/bash

flag=$1

[ "$flag" == "SKIP_INTEGRATION" ] && export SKIP_INTEGRATION=true 

[ "$SKIP_INTEGRATION" == true ] && echo "skiping integtation tests"

export PYTHONPATH=$PYTHONPATH:$PWD/src
python3 -m unittest
exit 0
