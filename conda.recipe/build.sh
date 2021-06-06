#!/bin/bash

set -euo pipefail

cd ${SRC_DIR}
IPOPT_DIR=${CONDA_PREFIX} ${PYTHON} -m pip install . -vv

cd ${SRC_DIR}/baseclasses
${PYTHON} -m pip install . -vv
