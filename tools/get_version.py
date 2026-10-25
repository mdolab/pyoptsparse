#!/usr/bin/env python3
# Run by meson.build to get the package version. This is the same as
# `python -m setuptools_scm --force-write-version-files`, but calling `python` by name would use whichever
# interpreter is first on PATH, which may not be the build environment's Python and may lack setuptools_scm.
#
# This file must stay NON-EXECUTABLE (git mode 100644). Meson only runs a non-executable script with the
# interpreter that is running meson (it rewrites the `/usr/bin/env python3` shebang). If the file is executable,
# the OS runs it via the shebang, `env` searches PATH for python3, and the original problem comes back.
# NumPy relies on the same behavior for numpy/_build_utils/gitversion.py.
import runpy
import sys

sys.argv = ["setuptools_scm", "--force-write-version-files"]
runpy.run_module("setuptools_scm", run_name="__main__")
