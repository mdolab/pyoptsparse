import os
import re
import shutil
import sys
import setuptools

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
    # check if ipopt dir is specified
    ipopt_dir_opt = ""
    if "IPOPT_DIR" in os.environ:
        ipopt_dir = os.environ["IPOPT_DIR"]
        ipopt_dir_opt = f"-Dipopt_dir={ipopt_dir}"
    prefix = os.path.join(os.getcwd(), staging_dir)
    purelibdir = "."

    # check if meson extra args are specified
    meson_args = ""
    if "MESON_ARGS" in os.environ:
        meson_args = os.environ["MESON_ARGS"]
    # check to make sure ipopt dir isnt specified twice
    if "-Dipopt_dir" in meson_args and ipopt_dir_opt != "":
        raise RuntimeError("IPOPT_DIR environment variable is set and '-Dipopt_dir' in MESON_ARGS")

    # configure
    meson_path = shutil.which("meson")
    meson_call = (
        f"{meson_path} setup {staging_dir} --prefix={prefix} "
        + f"-Dpython.purelibdir={purelibdir} -Dpython.platlibdir={purelibdir} {ipopt_dir_opt} {meson_args}"
    )
    sysargs = meson_call.split(" ")
    sysargs = [arg for arg in sysargs if arg != ""]
    meson_main(sysargs)

    # build
    meson_call = f"{meson_path} compile -C {staging_dir}"
    sysargs = meson_call.split(" ")
    meson_main(sysargs)


def copy_shared_libraries():
    build_path = os.path.join(staging_dir, "pyoptsparse")
    for root, _dirs, files in os.walk(build_path):
        for file in files:
            # move pyoptsparse to just under staging_dir
            if file.endswith((".so", ".lib", ".pyd", ".pdb", ".dylib")):
                if ".so.p" in root or ".pyd.p" in root:  # excludes intermediate object files
                    continue
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(staging_dir, new_path)
                new_path = new_path[match.span()[1] + 1:]
                print(f"Copying build file {file_path} -> {new_path}")
                shutil.copy(file_path, new_path)


if __name__ == "__main__":
    # This is where the meson build system will install to, it is then
    # used as the sources for setuptools
    staging_dir = "meson_build"

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
            "": ["*.so", "*.lib", "*.pyd", "*.pdb", "*.dylib", "assets/*", "LICENSE"],
        },
        python_requires=">=3.7",
        entry_points={
            "gui_scripts": [
                "optview = pyoptsparse.postprocessing.OptView:main",
                "optview_dash = pyoptsparse.postprocessing.OptView_dash:main",
            ]
        },
    )
