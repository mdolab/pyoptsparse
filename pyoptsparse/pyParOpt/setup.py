import os, sys


def configuration(parent_package="", top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration("pyParOpt", parent_package, top_path)

    return config
