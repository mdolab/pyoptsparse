#!/usr/bin/env python

import os,sys

def configuration(parent_package='',top_path=None):
	
    from numpy.distutils.misc_util import Configuration

    config = Configuration('pyNOMAD',parent_package,top_path)

    # platform-specific settings
    extra_link_args = []
    if sys.platform == "darwin":
	    extra_link_args.append("-bundle")
    #end
	
    config.add_library('nomad',
	    sources=[os.path.join('source/nomad_src', '*.cpp')])
    config.add_extension('_nomad',
	    sources=['source/NomadLinker.i','source/NomadLinker.cpp'],
	    include_dirs=['source', 'source/nomad_src'],
	    libraries=['nomad'],
	    extra_link_args=extra_link_args,
	    extra_compile_args=['-fPIC'],
	    swig_opts=[])
    config.add_data_files('LICENSE','README')

    return config
