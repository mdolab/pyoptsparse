#!/usr/local/bin/python

import os,sys

def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    # Check if we have source fils...if we don't we don't build:
    sources=[os.path.join('source', 'NLPQLP.F'),
             os.path.join('source', 'QL.F'),
             os.path.join('source', 'wrapper.F90')]
            
    config = Configuration('pyNLPQLP',parent_package,top_path)
    config.add_data_files('LICENSE','README')
    nlpqlp = os.path.join('pyoptsparse/pyNLPQLP/source', 'NLPQLP.F')

    if os.path.exists(nlpqlp):
        config.add_library('nlpqlp',sources=sources)
        config.add_extension('nlpqlp',sources=['source/f2py/nlpqlp.pyf'], libraries=['nlpqlp'])

    return config
    
