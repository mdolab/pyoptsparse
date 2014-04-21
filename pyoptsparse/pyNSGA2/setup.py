#!/usr/bin/env python

import os,sys

from numpy.distutils.command import build_src
from numpy.distutils.misc_util import appendpath
from numpy.distutils.command.build_src import get_swig_modulename, get_swig_target
from distutils.dep_util import newer_group
from numpy.distutils import log

#numpy/distutils/command/build_src.py
#584    if name != ext_name:
#587    ' but expected %r' % (source, name, ext_name))
#608    name = ext_name
#624    #py_files.append(os.path.join(py_target_dir, name+'.py'))
#640    swig_cmd = [swig, "-python", "-noproxy"]

def swig_sources(self, sources, extension):
	# Assuming SWIG 1.3.14 or later. See compatibility note in
	#   http://www.swig.org/Doc1.3/Python.html#Python_nn6
	
	new_sources = []
	swig_sources = []
	swig_targets = {}
	target_dirs = []
	py_files = []     # swig generated .py files
	target_ext = '.c'
	if self.swig_cpp:
		typ = 'c++'
		is_cpp = True
	else:
		typ = None
		is_cpp = False
	skip_swig = 0
	ext_name = extension.name.split('.')[-1]

	for source in sources:
		(base, ext) = os.path.splitext(source)
		if ext == '.i': # SWIG interface file
			if self.inplace:
				target_dir = os.path.dirname(base)
				py_target_dir = self.ext_target_dir
			else:
				target_dir = appendpath(self.build_src, os.path.dirname(base))
				py_target_dir = target_dir
			if os.path.isfile(source):
				name = get_swig_modulename(source)
				if name != ext_name:
					raise DistutilsSetupError(
						'mismatch of extension names: %s provides %r'
						' but expected %r' % (source, name, ext_name))
				if typ is None:
					typ = get_swig_target(source)
					is_cpp = typ=='c++'
					if is_cpp: target_ext = '.cpp'
				else:
					typ2 = get_swig_target(source)
					if typ!=typ2:
						log.warn('expected %r but source %r defines %r swig target' \
								% (typ, source, typ2))
						if typ2=='c++':
							log.warn('resetting swig target to c++ (some targets may have .c extension)')
							is_cpp = True
							target_ext = '.cpp'
						else:
							log.warn('assuming that %r has c++ swig target' % (source))
				target_file = os.path.join(target_dir,'%s_wrap%s' \
										% (name, target_ext))
			else:
				log.warn('  source %s does not exist: skipping swig\'ing.' \
						% (source))
				name = ext_name
				skip_swig = 1
				target_file = _find_swig_target(target_dir, name)
				if not os.path.isfile(target_file):
					log.warn('  target %s does not exist:\n   '\
							'Assuming %s_wrap.{c,cpp} was generated with '\
							'"build_src --inplace" command.' \
							% (target_file, name))
					target_dir = os.path.dirname(base)
					target_file = _find_swig_target(target_dir, name)
					if not os.path.isfile(target_file):
						raise DistutilsSetupError("%r missing" % (target_file,))
					log.warn('   Yes! Using %r as up-to-date target.' \
							% (target_file))
			target_dirs.append(target_dir)
			new_sources.append(target_file)
			#py_files.append(os.path.join(py_target_dir, name+'.py'))
			swig_sources.append(source)
			swig_targets[source] = new_sources[-1]
		else:
			new_sources.append(source)

	if not swig_sources:
		return new_sources

	if skip_swig:
		return new_sources + py_files

	for d in target_dirs:
		self.mkpath(d)

	swig = self.swig or self.find_swig()
	swig_cmd = [swig, "-python", "-noproxy"]
	if is_cpp:
		swig_cmd.append('-c++')
	for d in extension.include_dirs:
		swig_cmd.append('-I'+d)
	for source in swig_sources:
		target = swig_targets[source]
		depends = [source] + extension.depends
		if self.force or newer_group(depends, target, 'newer'):
			log.info("%s: %s" % (os.path.basename(swig) \
								+ (is_cpp and '++' or ''), source))
			self.spawn(swig_cmd + self.swig_opts \
					+ ["-o", target, '-outdir', py_target_dir, source])
		else:
			log.debug("  skipping '%s' swig interface (up-to-date)" \
					% (source))

	return new_sources + py_files
	

build_src.build_src.swig_sources = swig_sources


def configuration(parent_package='',top_path=None):
	
	from numpy.distutils.misc_util import Configuration
	
	config = Configuration('pyNSGA2',parent_package,top_path)
	
	# platform-specific settings
	extra_link_args = []
	if sys.platform == "darwin":
		extra_link_args.append("-bundle")
	#end
	
	config.add_library('nsga2',
		sources=[os.path.join('source', '*.c')])
	config.add_extension('nsga2',
		sources=['source/swig/nsga2.i'],
		include_dirs=['source'],
		libraries=['nsga2'],
		extra_link_args=extra_link_args)
	config.add_data_files('LICENSE','README')
	
	return config
	

if __name__ == '__main__':
	from numpy.distutils.core import setup
	setup(**configuration(top_path='').todict())
	
