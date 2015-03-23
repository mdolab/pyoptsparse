#!/usr/local/bin/python

import os,sys


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyFSQP',parent_package,top_path)
    
    config.add_library('ffsqp',
        sources=[os.path.join('source', '*.f')])
    config.add_extension('ffsqp',
        sources=['source/f2py/ffsqp.pyf'],
        libraries=['ffsqp'])
    config.add_data_files('LICENSE','README')
    
    return config
