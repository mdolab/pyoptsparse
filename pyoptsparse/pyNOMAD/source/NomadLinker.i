//-*- c++ -*-
%module nomad
%module(directors="1") nomad 
%{
    #define SWIG_FILE_WITH_INIT
    #include "Callback.h"
    #include "NomadLinker.h"
%}
%feature("director") Callback;
%feature("nodirector") NomadLinker;

%include "std_vector.i"
namespace std {
    %template(fvector) vector<double>;
}
    
%include "numpy.i"
%init %{
import_array();
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* x, int xdim)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* lb, int lbdim)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double* ub, int ubdim)}

%feature("pythonprepend") NomadLinker::setCallback(Callback&) %{
   if len(args) == 1 and (not isinstance(args[0], Callback) and callable(args[0])):
      class CallableWrapper(Callback):
         def __init__(self, f):
            super(CallableWrapper, self).__init__()
            self.f_ = f
         def call(self, obj, x):
            fx = self.f_(obj, x)
            obj.setSolution(fx)
      args = tuple([CallableWrapper(args[0])])
      args[0].__disown__()
   elif len(args) == 1 and isinstance(args[0], Callback):
      args[0].__disown__()


%}

%include "Callback.h"
%include "NomadLinker.h"
