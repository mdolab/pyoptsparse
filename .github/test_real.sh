#!/bin/bash
set -e

# all tests should pass on the private image
# except for the Intel image, where IPOPT is not available
if [[ $IMAGE == "private" ]]; then
    EXTRA_FLAGS='--disallow_skipped'
fi

# Set OpenMPI env variables only on non-Intel MPI
if [[ -z $I_MPI_ROOT ]]; then
    # Set these to allow MPI oversubscription because the tests need to run on specific number of procs but the test runner may have fewer
    export OMPI_MCA_rmaps_base_oversubscribe=1 # This works for OpenMPI <= 4
    export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe # This works from OpenMPI >= 5
fi

cd tests
# we have to copy over the coveragerc file to make sure it's in the
# same directory where codecov is run
cp ../.coveragerc .
testflo -i --pre_announce --disallow_deprecations -v --coverage --coverpkg pyoptsparse $EXTRA_FLAGS --timeout 120
