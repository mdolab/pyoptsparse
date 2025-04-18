import os
import re
import shutil
import setuptools
import subprocess


def run_meson_build():
    prefix = os.path.join(os.getcwd(), staging_dir)
    purelibdir = "."

    # check if meson extra args are specified
    meson_args = ""
    if "MESON_ARGS" in os.environ:
        meson_args = os.environ["MESON_ARGS"]

    # configure
    meson_path = shutil.which("meson")
    meson_call = (
        f"{meson_path} setup {staging_dir} --prefix={prefix} "
        + f"-Dpython.purelibdir={purelibdir} -Dpython.platlibdir={purelibdir} {meson_args}"
    )
    sysargs = meson_call.split(" ")
    sysargs = [arg for arg in sysargs if arg != ""]
    p1 = subprocess.run(sysargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p1.stdout.decode())
    setup_log = os.path.join(staging_dir, "setup.log")
    with open(setup_log, "wb") as f:
        f.write(p1.stdout)
    if p1.returncode != 0:
        raise OSError(sysargs, f"The meson setup command failed! Check the log at {setup_log} for more information.")

    # build
    meson_call = f"{meson_path} compile -C {staging_dir}"
    sysargs = meson_call.split(" ")
    p2 = subprocess.run(sysargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p2.stdout.decode())
    compile_log = os.path.join(staging_dir, "compile.log")
    with open(compile_log, "wb") as f:
        f.write(p2.stdout)
    if p2.returncode != 0:
        raise OSError(
            sysargs, f"The meson compile command failed! Check the log at {compile_log} for more information."
        )


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
                new_path = new_path[match.span()[1] + 1 :]
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
            "packaging",
            "sqlitedict>=1.6",
            "numpy>=1.21",
            "scipy>=1.7",
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
        python_requires=">=3.9",
        entry_points={
            "gui_scripts": [
                "optview = pyoptsparse.postprocessing.OptView:main",
                "optview_dash = pyoptsparse.postprocessing.OptView_dash:main",
            ]
        },
    )
