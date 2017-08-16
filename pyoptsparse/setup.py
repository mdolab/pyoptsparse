#!/usr/bin/env python

import os,sys

def configuration(parent_package='', top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyoptsparse', parent_package, top_path)
    
    # need: auto add_subpackage from source availability
    config.add_subpackage('pySNOPT')
    config.add_subpackage('pyIPOPT')
    config.add_subpackage('pySLSQP')
    config.add_subpackage('pyCONMIN')
    config.add_subpackage('pyFSQP')
    config.add_subpackage('pyNLPQLP')
    config.add_subpackage('pyNSGA2')
    config.add_subpackage('pyPSQP')
    config.add_subpackage('pyNLPY_AUGLAG')
    config.add_subpackage('pyALPSO')
#    config.add_subpackage('pyNOMAD')
    config.add_subpackage('sqlitedict')
    config.add_data_files('LICENSE')
    
    return config
 
