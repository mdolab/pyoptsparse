"""Compute the package version with setuptools_scm and write pyoptsparse/_version.py.

Called by meson at configure time (see project() in meson.build). The version comes from,
in order of priority: the SETUPTOOLS_SCM_PRETEND_VERSION environment variable (set e.g. by
the conda-forge recipe when building from a GitHub tarball), git tags, or PKG-INFO when
building from an sdist.
"""

# External modules
from setuptools_scm import get_version

if __name__ == "__main__":
    # root/version_file must be relative paths, anchored to this file's directory via relative_to
    print(get_version(root=".", relative_to=__file__, version_file="pyoptsparse/_version.py"))
