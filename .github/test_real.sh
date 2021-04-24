#!/bin/bash
set -e
if [[ $IMAGE == "public" ]]; then
    testflo --pre_announce -v
else
    # all tests should pass on private
    testflo --pre_announce -v --disallow_skipped
fi
