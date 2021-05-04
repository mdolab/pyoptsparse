#!/bin/bash
set -e

# all tests should pass on private
if [[ $IMAGE == "private" ]]; then
    EXTRA_FLAGS='--disallow_skipped'
fi

testflo --pre_announce -v --coverage --coverpkg pyoptsparse $EXTRA_FLAGS
