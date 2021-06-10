#!/bin/bash
set -e

# first copy back proprietary solvers
if [[ $IMAGE == "private" ]]; then
    cp -r $HOME/NLPQLP/* pyoptsparse/pyNLPQLP/source
    cp -r $HOME/SNOPT/* pyoptsparse/pySNOPT/source
fi

# install in place due to coverage
pip install -e .[optview,testing]
