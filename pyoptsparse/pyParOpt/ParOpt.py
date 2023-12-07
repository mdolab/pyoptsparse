try:
    from paropt.paropt_pyoptsparse import ParOptSparse as ParOpt
except:

    class ParOpt:
        def __init__(self, raiseError=True, options={}):
            name = "ParOpt"
            category = "Local Optimizer"
            if raiseError:
                raise Error("There was an error importing ParOpt")
