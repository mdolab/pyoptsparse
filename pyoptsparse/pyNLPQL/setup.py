#!/usr/local/bin/python

import os,sys


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyNLPQL',parent_package,top_path)
    
    config.add_library('nlpql',
        sources=[os.path.join('source', '*.f')])
    config.add_extension('nlpql',
        sources=['source/f2py/nlpql.pyf'],
        libraries=['nlpql'])
    config.add_data_files('LICENSE','README')
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
    
