#!/usr/local/bin/python

import os,sys
from numpy.distutils.command.build_ext import build_ext

if os.path.exists('MANIFEST'): 
    os.remove('MANIFEST')
#end

if sys.version_info[:2] < (2, 5):
    print('pyOptSparse requires Python version 2.5 or later (%d.%d detected).' %sys.version_info[:2])
    sys.exit(-1)
#end

try:                  
    import numpy
    if int(numpy.__version__.split('.')[0]) < 1:
        print('pyOptSparse requires NumPy version 1.0 or later (%s detected).' %numpy.__version__)
        sys.exit(-1)
    #end
except ImportError as e:
    print('NumPy version 1.0 or later must be installed to build pyOptSparse')
    sys.exit(-1)
#end

if sys.argv[-1].endswith('setup.py'):
    print('\nTo install, run "python setup.py install"\n\nTo build, run "python setup.py inplace"\n')
    sys.exit(-1)
#end

if sys.argv[1] == 'inplace':
    sys.argv[1:2] = ['build_src','--inplace','build_ext','--inplace','build']
#end

if sys.argv[1] == 'install':
    arg = []
    i = 2
    while i < len(sys.argv):
        if (sys.argv[i].startswith('--compiler=') or sys.argv[i].startswith('--fcompiler=')):
            arg.append(sys.argv[i])
            del sys.argv[i]
        else:
            i += 1
        #end
    #end
    if arg != []:
        arg.insert(0,'build')
        sys.argv[1:1] = arg
    #end
#end

if sys.argv[1] == 'compilers':
    del sys.argv[1]
    from numpy.distutils.fcompiler import show_fcompilers
    show_fcompilers()
    print('\n')
    from numpy.distutils.ccompiler import show_compilers
    show_compilers()
    sys.exit()
#end


class build_opt(build_ext):
    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except:
            self.announce('*** WARNING: Building of optimizer "%s" '
            'failed: %s' %(ext.name, sys.exc_info()[1]))
        #end
    #end
#end


def configuration(parent_package='',top_path=None):
    
    from numpy.distutils.misc_util import Configuration
    
    config = Configuration(None,parent_package,top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    
    config.add_subpackage('pyOptSparse')
    
    return config
    

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        name             = 'pyOptSparse',
        version          = '1.0.0',
        author           = 'Dr. Gaetan Kenway',
        author_email     = 'gaetank@gmail.com',
        maintainer       = 'Dr. Gaetan Kenway',
        maintainer_email = 'gaetank@gmail.com',
        url              = 'http://pyopt.org/',
        download_url     = 'http://pyopt.org/',
        description      = 'Python package for formulating and solving nonlinear constrained optimization problems',
        long_description = 'pyOpt is a Python package for formulating and solving nonlinear constrained optimization problems',
        keywords         = 'optimization',
        license          = 'GNU LGPL',
        platforms        = ['Windows','Linux','Solaris','Mac OS-X','Unix'],
        classifiers      = ['Development Status :: 5 - Production/Stable',
                            'Environment :: Console',
                            'Intended Audience :: Science/Research',
                            'Intended Audience :: Developers',
                            'Intended Audience :: Education',
                            'License :: LGPL',
                            'Operating System :: Microsoft :: Windows',
                            'Operating System :: POSIX :: Linux',
                            'Operating System :: Unix',
                            'Operating System :: MacOS',
                            'Programming Language :: Python',
                            'Topic :: Scientific/Engineering',
                            'Topic :: Software Development',
                            'Topic :: Education'],
        configuration    = configuration,
        cmdclass = {"build_ext": build_opt},
        #**configuration().todict()
    )
    
