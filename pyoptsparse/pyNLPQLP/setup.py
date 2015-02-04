#!/usr/local/bin/python

import os,sys


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration('pyNLPQLP',parent_package,top_path)
    
    config.add_library('nlpqlp',
        sources=[os.path.join('source', 'NLPQLP.F'),
                 os.path.join('source', 'QL.F'),
                 os.path.join('source', 'wrapper.F90')])
    config.add_extension('nlpqlp',
        sources=['source/f2py/nlpqlp.pyf'],
        libraries=['nlpqlp'])
    config.add_data_files('LICENSE','README')
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
