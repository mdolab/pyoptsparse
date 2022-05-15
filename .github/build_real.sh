#!/bin/bash
set -e

# first copy back proprietary solvers
if [[ $IMAGE == "private" ]]; then
    cp -r $HOME/NLPQLP/* pyoptsparse/pyNLPQLP/source
    cp -r $HOME/SNOPT/* pyoptsparse/pySNOPT/source
fi

export IPOPT_DIR="pyoptsparse/pyIPOPT/Ipopt"
pip install .[optview,testing]
