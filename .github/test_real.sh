#!/bin/bash
set -e

# all tests should pass on the private image
# except for the Intel image, where IPOPT is not available
if [[ $IMAGE == "private" ]]; then
    EXTRA_FLAGS='--disallow_skipped'
fi

cd tests
# we have to copy over the coveragerc file to make sure it's in the
# same directory where codecov is run
cp ../.coveragerc .
testflo -i --pre_announce --disallow_deprecations -v --coverage --coverpkg pyoptsparse $EXTRA_FLAGS
