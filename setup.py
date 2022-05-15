import os
import re
import shutil
import sys

import setuptools  # magic import to allow us to use entry_point


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
        ipopt_dir_opt = f'-Dipopt_dir="{ipopt_dir}"'
    os.system(f"meson setup meson_build --prefix" f"={os.path.join(os.getcwd(), 'staging_dir')}" f" {ipopt_dir_opt}")
    os.system(f"ninja -C meson_build install")
    use_local_install_dir()
    copy_shared_libraries()


def use_local_install_dir():
    installation = ""
    for root, dirs, files in os.walk("staging_dir"):
        for dir in dirs:

            # move pyoptsparse to just under staging_dir
            if dir.endswith("pyoptsparse"):
                installation = os.path.join(root, dir)
    new_installation = os.path.join("staging_dir", "pyoptsparse")
    shutil.move(installation, new_installation)
    shutil.rmtree(os.path.join("staging_dir", "lib"))


def copy_shared_libraries():
    so_files = []
    for root, dirs, files in os.walk("staging_dir"):
        for file in files:

            # move pyoptsparse to just under staging_dir
            if file.endswith(".so"):
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(staging_dir, new_path)
                new_path = new_path[match.span()[1] + 1 :]
                shutil.copy(file_path, new_path)


if __name__ == "__main__":
    # This is where the meson build system will install to, it is then
    # used as the sources for setuptools
    staging_dir = "staging_dir"

    doc_path = os.path.join(staging_dir, "doc")

    # this keeps the meson build system from running more than once
    if "dist" not in str(os.path.abspath(__file__)) and not os.path.isdir(staging_dir):
        run_meson_build()
        shutil.copytree("doc", doc_path)

    req_txt = os.path.join(doc_path, "requirements.txt")
    with open(req_txt) as f:
        print(f"THIS IS REQ PATH: {req_txt}")
        docs_require = f.read().splitlines()

    init_file = os.path.join(staging_dir, "pyoptsparse", "__init__.py")
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
            "testing": ["testflo>=1.4.5"],
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
        package_dir={"pyoptsparse": os.path.join(staging_dir, "pyoptsparse"), "doc": os.path.join(staging_dir, "doc")},
        # packages=setuptools.find_packages(where=staging_dir),
        python_requires=">=3.9",
        package_data={
            "pyoptsparse": ["*.so", "*.lib", "assets/*"],
            "doc": ["*"],
        },
        entry_points={
            "gui_scripts": [
                "optview = pyoptsparse.postprocessing.OptView:main",
                "optview_dash = pyoptsparse.postprocessing.OptView_dash:main",
            ]
        },
    )
