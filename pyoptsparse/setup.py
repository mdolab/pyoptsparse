#!/usr/bin/env python

import os,sys

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyoptsparse',parent_package,top_path)
    
    # need: auto add_subpackage from source availability
    config.add_subpackage('pySNOPT')
    config.add_subpackage('pyIPOPT')
    config.add_subpackage('pySLSQP')
    config.add_subpackage('pyCONMIN')
    config.add_subpackage('pyFSQP')
    config.add_subpackage('pyNLPQL')
    config.add_subpackage('pyNSGA2')
    config.add_subpackage('pyPSQP')
    config.add_subpackage('pyNLPY_AUGLAG')
    config.add_data_files('LICENSE','README')
    
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
    
