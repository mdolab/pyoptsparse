#!/bin/bash
set -e

# all tests should pass on the private image
if [[ $IMAGE == "private" ]]; then
    EXTRA_FLAGS='--disallow_skipped'
fi

cd tests
cp ../.coveragerc .
testflo --pre_announce -v --coverage --coverpkg pyoptsparse $EXTRA_FLAGS
