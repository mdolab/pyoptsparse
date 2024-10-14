#!/bin/bash
set -e

# first copy back proprietary solvers
if [[ $IMAGE == "private" ]]; then
    cp -r $HOME/NLPQLP/* pyoptsparse/pyNLPQLP/source
    cp -r $HOME/SNOPT/* pyoptsparse/pySNOPT/source
fi

# temporarily disable pip constraints file for build-time dependencies
unset PIP_CONSTRAINT

pip install .[optview,testing] -v
