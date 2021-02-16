# back up proprietary solvers
if [[ $IMAGE == "private" ]]; then
    cp -r pyoptsparse/pySNOPT/source $HOME/SNOPT
    cp -r pyoptsparse/pyNLPQLP/source $HOME/NLPQLP
fi
