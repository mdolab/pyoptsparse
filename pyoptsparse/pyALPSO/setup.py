import os, sys


def configuration(parent_package="", top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration("pyALPSO", parent_package, top_path)

    config.add_data_files("LICENSE", "README")

    return config
