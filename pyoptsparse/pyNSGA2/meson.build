swig = find_program('swig', required: true)

nsga2_source = custom_target('nsga2_wrap.c',
                            input : ['source/swig/nsga2.i'],
                            output : ['nsga2_wrap.c'],
                            command: ['swig', '-o', 'pyoptsparse/pyNSGA2/nsga2_wrap.c', '-python', '-interface', 'nsga2', '@INPUT@']
                            )

py3_target.extension_module('nsga2',
                'source/allocate.c',
                'source/auxiliary.c',
                'source/crossover.c',
                'source/crowddist.c',
                'source/decode.c',
                'source/dominance.c',
                'source/eval.c',
                'source/fillnds.c',
                'source/initialize.c',
                'source/list.c',
                'source/merge.c',
                'source/mutation.c',
                'source/nsga2.c',
                'source/rand.c',
                'source/rank.c',
                'source/report.c',
                'source/sort.c',
                'source/tourselect.c',
                nsga2_source,
                include_directories: 'source',
                dependencies : py3_dep,
                subdir: 'pyoptsparse/pyNSGA2',
                link_language: 'c',
                install : false)

#python_sources = [
#    '__init__.py',
#    'pyNSGA2.py',
#    'LICENSE'
#]
#
#py3_target.install_sources(
#  python_sources,
#  pure: false,
#  subdir: 'pyoptsparse/pyNSGA2'
#)
