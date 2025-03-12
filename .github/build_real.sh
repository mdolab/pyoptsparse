#!/bin/bash
set -e

# first copy back proprietary solvers
if [[ $IMAGE == "private" ]]; then
    cp -r $HOME/NLPQLP/* pyoptsparse/pyNLPQLP/source
    cp -r $HOME/SNOPT/* pyoptsparse/pySNOPT/source
fi

# temporarily disable pip constraints file for build-time dependencies
mv ~/.config/pip/constraints.txt ~/.config/pip/constraints.txt.bkup
touch  ~/.config/pip/constraints.txt

# set $PKG_CONFIG_PATH so pkg-config can find ipopt
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$IPOPT_DIR/lib/pkgconfig

pip install .[optview,testing] -v

# move pip constraints file back
mv ~/.config/pip/constraints.txt.bkup ~/.config/pip/constraints.txt