#!/bin/bash

set -euo pipefail

cd ${SRC_DIR}
${PYTHON} -m pip install . -vv

cd ${SRC_DIR}/baseclasses
${PYTHON} -m pip install . -vv
