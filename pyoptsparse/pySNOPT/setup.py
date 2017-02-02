#!/usr/bin/env python

import os,sys


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pySNOPT', parent_package, top_path)
    config.add_data_files('LICENSE','README')

    # Since snopt has a bunch of source files, we will just check if
    # snoptc.c exists. If so, we will assume all the rest of the files
    # are present. 
    
    snoptc = os.path.join('pyoptsparse/pySNOPT/source', 'snoptc.f')
    if os.path.exists(snoptc):
        config.add_library('snopt', sources=[os.path.join('source', '*.f')])
        config.add_extension('snopt', sources=['source/f2py/snopt.pyf'],
                             libraries=['snopt'])
    return config
   
