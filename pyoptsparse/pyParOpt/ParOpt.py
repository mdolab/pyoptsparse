try:
    from paropt.paropt_pyoptsparse import ParOptSparse as ParOpt
except:

    class ParOpt:
        pass
