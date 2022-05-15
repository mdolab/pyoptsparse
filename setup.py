import copy
import os
import re
import shutil
import sys

import setuptools  # magic import to allow us to use entry_point

from mesonbuild.mesonmain import run as meson_run


def meson_main(sysargs):
    # Always resolve the command path so Ninja can find it for regen, tests, etc.
    if "meson.exe" in sys.executable:
        assert os.path.isabs(sys.executable)
        launcher = sys.executable
    else:
        launcher = os.path.realpath(sysargs[0])
    return meson_run(sysargs[1:], launcher)


def run_meson_build():
    # set compilers for specific platforms
    if sys.platform == "win32":
        os.environ["CC"] = "cl"
        os.environ["FC"] = "flang"
        os.environ["CC_LD"] = "link"
    elif sys.platform == "darwin":
        os.environ["CC"] = "clang"
        os.environ["FC"] = "gfortran"
        os.environ["CC_LD"] = "ld"
    # elif sys.platform == "linux":
    # os.environ["CC"] = "gcc"
    # os.environ["FC"] = "gfortran"
    # os.environ["CC_LD"] = "bfd"

    # check if ipopt dir is specified
    ipopt_dir_opt = ""
    if "IPOPT_DIR" in os.environ:
        ipopt_dir = os.environ["IPOPT_DIR"]
        ipopt_dir_opt = f'-Dipopt_dir={ipopt_dir}'
    prefix = os.path.join(os.getcwd(), staging_dir)
    purelibdir = "."

    # configure
    meson_path = shutil.which("meson")
    meson_call = (
        f"{meson_path} setup meson_build --prefix={prefix} "
        f"-Dpython.purelibdir={purelibdir} -Dpython.platlibdir={purelibdir} {ipopt_dir_opt}"
    )
    sysargs = meson_call.split(" ")
    sysargs = [arg for arg in sysargs if arg != ""]
    meson_main(sysargs)

    # build
    meson_call = f"{meson_path} install -C meson_build"
    sysargs = meson_call.split(" ")
    meson_main(sysargs)


def copy_shared_libraries():
    for root, dirs, files in os.walk(staging_dir):
        for file in files:

            # move pyoptsparse to just under staging_dir
            if file.endswith((".so", ".lib", ".pyd", ".pdb", ".dylib")):
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(staging_dir, new_path)
                new_path = new_path[match.span()[1] + 1 :]
                shutil.copy(file_path, new_path)


if __name__ == "__main__":
    # This is where the meson build system will install to, it is then
    # used as the sources for setuptools
    staging_dir = "staging_dir"

    # this keeps the meson build system from running more than once
    if "dist" not in str(os.path.abspath(__file__)) and not os.path.isdir(staging_dir):
        cwd = os.getcwd()
        run_meson_build()
        os.chdir(cwd)
        copy_shared_libraries()

    docs_require = ""
    req_txt = os.path.join("doc", "requirements.txt")
    if os.path.isfile(req_txt):
        with open(req_txt) as f:
            docs_require = f.read().splitlines()

    init_file = os.path.join("pyoptsparse", "__init__.py")
    __version__ = re.findall(
        r"""__version__ = ["']+([0-9\.]*)["']+""",
        open(init_file).read(),
    )[0]

    setuptools.setup(
        name="pyoptsparse",
        version=__version__,
        description="Python package for formulating and solving nonlinear constrained optimization problems",
        long_description="pyOptSparse is a Python package for formulating and solving nonlinear constrained optimization problems",
        platforms=["Linux"],
        keywords="optimization",
        install_requires=[
            "sqlitedict>=1.6",
            "numpy>=1.16",
            "scipy>1.2",
            "mdolab-baseclasses>=1.3.1",
        ],
        extras_require={
            "optview": [
                "dash",
                "plotly",
                "matplotlib",
            ],
            "docs": docs_require,
            "testing": ["testflo>=1.4.5", "parameterized"],
        },
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "Intended Audience :: Education",
            "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering",
            "Topic :: Software Development",
            "Topic :: Education",
        ],
        package_dir={"": "."},
        packages=setuptools.find_packages(where="."),
        package_data={
            "": ["*.so", "*.lib", "*.pyd", "*.pdb", "*.dylib", "assets/*"],
        },
        python_requires=">=3.8",
        entry_points={
            "gui_scripts": [
                "optview = pyoptsparse.postprocessing.OptView:main",
                "optview_dash = pyoptsparse.postprocessing.OptView_dash:main",
            ]
        },
    )
