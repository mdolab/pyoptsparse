/*
 * Copyright (c) 2008, Eric You Xu, Washington University All rights
 * reserved. Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following conditions
 * are met:
 * 
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. * Redistributions in
 * binary form must reproduce the above copyright notice, this list of
 * conditions and the following disclaimer in the documentation and/or other
 * materials provided with the distribution. * Neither the name of the
 * Washington University nor the names of its contributors may be used to
 * endorse or promote products derived from this software without specific
 * prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

/* 
 * Added "eval_intermediate_callback" by 
 * OpenMDAO at NASA Glenn Research Center, 2010 and 2011
 *
 * Changed logger from code contributed by alanfalloon  
*/

#include "hook.h"
#include <unistd.h>

void logger(const char *fmt, ...)
{
	if (user_log_level == VERBOSE) {
		va_list ap;
		va_start(ap, fmt);
		PySys_WriteStdout(fmt, ap);
		va_end(ap);
		PySys_WriteStdout("\n");
	}
}

Bool eval_intermediate_callback(Index alg_mod,	/* 0 is regular, 1 is resto */
				Index iter_count, Number obj_value,
				Number inf_pr, Number inf_du,
				Number mu, Number d_norm,
				Number regularization_size,
				Number alpha_du, Number alpha_pr,
				Index ls_trials, UserDataPtr data)
{
	//logger("[Callback:E]intermediate_callback");

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	long result_as_long;
	Bool result_as_bool;

	PyObject *python_algmod = Py_BuildValue("i", alg_mod);
	PyObject *python_iter_count = Py_BuildValue("i", iter_count);
	PyObject *python_obj_value = Py_BuildValue("d", obj_value);
	PyObject *python_inf_pr = Py_BuildValue("d", inf_pr);
	PyObject *python_inf_du = Py_BuildValue("d", inf_du);
	PyObject *python_mu = Py_BuildValue("d", mu);
	PyObject *python_d_norm = Py_BuildValue("d", d_norm);
	PyObject *python_regularization_size =
	    Py_BuildValue("d", regularization_size);
	PyObject *python_alpha_du = Py_BuildValue("d", alpha_du);
	PyObject *python_alpha_pr = Py_BuildValue("d", alpha_pr);
	PyObject *python_ls_trials = Py_BuildValue("i", ls_trials);

	PyObject *arglist = NULL;

	if (user_data != NULL)
		arglist = Py_BuildValue("(OOOOOOOOOOOO)",
					python_algmod,
					python_iter_count,
					python_obj_value,
					python_inf_pr,
					python_inf_du,
					python_mu,
					python_d_norm,
					python_regularization_size,
					python_alpha_du,
					python_alpha_pr,
					python_ls_trials,
					(PyObject *) user_data);
	else
		arglist = Py_BuildValue("(OOOOOOOOOOO)",
					python_algmod,
					python_iter_count,
					python_obj_value,
					python_inf_pr,
					python_inf_du,
					python_mu,
					python_d_norm,
					python_regularization_size,
					python_alpha_du,
					python_alpha_pr, python_ls_trials);

	PyObject *result =
	    PyObject_CallObject(myowndata->eval_intermediate_callback_python,
				arglist);

	if (!result)
		PyErr_Print();

	result_as_long = PyLong_AsLong(result);
	result_as_bool = (Bool) result_as_long;

	Py_DECREF(result);
	Py_CLEAR(arglist);
	//logger("[Callback:R] intermediate_callback");
	return result_as_bool;
}

Bool
eval_f(Index n, Number * x, Bool new_x, Number * obj_value, UserDataPtr data)
{
	//logger("[Callback:E] eval_f");

	npy_intp dims[1];
	dims[0] = n;

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	// import_array ();

	import_array1(FALSE);
	PyObject *arrayx =
	    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE, (char *)x);
	if (!arrayx)
		return FALSE;

	if (new_x && myowndata->apply_new_python) {
		/* Call the python function to applynew */
		PyObject *arg1;
		arg1 = Py_BuildValue("(O)", arrayx);
		PyObject *tempresult = PyObject_CallObject(
        myowndata->apply_new_python, arg1);
		if (tempresult == NULL) {
			logger("[Error] Python function apply_new returns NULL");
      PyErr_Print();
			Py_DECREF(arg1);
			return FALSE;
		}
		Py_DECREF(arg1);
		Py_DECREF(tempresult);
	}

	PyObject *arglist;
	if (user_data != NULL) {
		arglist = Py_BuildValue("(OO)", arrayx, (PyObject *) user_data);
  } else {
		arglist = Py_BuildValue("(O)", arrayx);
  }

	PyObject *result = PyObject_CallObject(myowndata->eval_f_python, arglist);

	if (result == NULL) {
    logger("[Error] Python function eval_f returns NULL");
		PyErr_Print();
		Py_DECREF(arrayx);
		Py_CLEAR(arglist);
		return FALSE;
	}

	*obj_value = PyFloat_AsDouble(result);

  if (PyErr_Occurred()) {
    logger("[Error] Python function eval_f returns non-PyFloat");
		PyErr_Print();
		Py_DECREF(result);
		Py_DECREF(arrayx);
		Py_CLEAR(arglist);
		return FALSE;
  }

	Py_DECREF(result);
	Py_DECREF(arrayx);
	Py_CLEAR(arglist);
	//logger("[Callback:R] eval_f");
	return TRUE;
}

Bool
eval_grad_f(Index n, Number * x, Bool new_x, Number * grad_f, UserDataPtr data)
{
	//logger("[Callback:E] eval_grad_f");

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	if (myowndata->eval_grad_f_python == NULL)
		PyErr_Print();

	/* int dims[1]; */
	npy_intp dims[1];
	dims[0] = n;
	// import_array ();

	import_array1(FALSE);

	/*
	 * PyObject *arrayx = PyArray_FromDimsAndData(1, dims, PyArray_DOUBLE
	 * , (char*) x);
	 */
	PyObject *arrayx =
	    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE, (char *)x);
	if (!arrayx)
		return FALSE;

	if (new_x && myowndata->apply_new_python) {
		/* Call the python function to applynew */
		PyObject *arg1 = Py_BuildValue("(O)", arrayx);
		PyObject *tempresult = PyObject_CallObject(
        myowndata->apply_new_python, arg1);
		if (tempresult == NULL) {
			logger("[Error] Python function apply_new returns NULL");
      PyErr_Print();
			Py_DECREF(arg1);
			return FALSE;
		}
		Py_DECREF(arg1);
		Py_DECREF(tempresult);
	}

	PyObject *arglist;
	if (user_data != NULL)
		arglist = Py_BuildValue("(OO)", arrayx, (PyObject *) user_data);
	else
		arglist = Py_BuildValue("(O)", arrayx);

	PyArrayObject *result = (PyArrayObject *) PyObject_CallObject(
      myowndata->eval_grad_f_python, arglist);

	if (result == NULL) {
    logger("[Error] Python function eval_grad_f returns NULL");
		PyErr_Print();
    return FALSE;
  }
  
  if (!PyArray_Check(result)) {
    logger("[Error] Python function eval_grad_f returns non-PyArray");
    Py_DECREF(result);
    return FALSE;
  }

	double *tempdata = (double *)result->data;
	int i;
	for (i = 0; i < n; i++)
		grad_f[i] = tempdata[i];

	Py_DECREF(result);
	Py_CLEAR(arrayx);
	Py_CLEAR(arglist);
	//logger("[Callback:R] eval_grad_f");
	return TRUE;
}

Bool
eval_g(Index n, Number * x, Bool new_x, Index m, Number * g, UserDataPtr data)
{

	//logger("[Callback:E] eval_g");

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	if (myowndata->eval_g_python == NULL)
		PyErr_Print();
	/* int dims[1]; */
	npy_intp dims[1];
	int i;
	double *tempdata;

	dims[0] = n;
	// import_array ();

	import_array1(FALSE);

	/*
	 * PyObject *arrayx = PyArray_FromDimsAndData(1, dims, PyArray_DOUBLE
	 * , (char*) x);
	 */
	PyObject *arrayx =
	    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE, (char *)x);
	if (!arrayx)
		return FALSE;

	if (new_x && myowndata->apply_new_python) {
		/* Call the python function to applynew */
		PyObject *arg1 = Py_BuildValue("(O)", arrayx);
		PyObject *tempresult = PyObject_CallObject(
        myowndata->apply_new_python, arg1);
		if (tempresult == NULL) {
			logger("[Error] Python function apply_new returns NULL");
      PyErr_Print();
			Py_DECREF(arg1);
			return FALSE;
		}
		Py_DECREF(arg1);
		Py_DECREF(tempresult);
	}

	PyObject *arglist;
	if (user_data != NULL)
		arglist = Py_BuildValue("(OO)", arrayx, (PyObject *) user_data);
	else
		arglist = Py_BuildValue("(O)", arrayx);

	PyArrayObject *result = (PyArrayObject *) PyObject_CallObject(
      myowndata->eval_g_python, arglist);

  if (result == NULL) {
    logger("[Error] Python function eval_g returns NULL");
		PyErr_Print();
    return FALSE;
  }
  
  if (!PyArray_Check(result)) {
    logger("[Error] Python function eval_g returns non-PyArray");
    Py_DECREF(result);
    return FALSE;
  }

	tempdata = (double *)result->data;
	for (i = 0; i < m; i++) {
		g[i] = tempdata[i];
	}

	Py_DECREF(result);
	Py_CLEAR(arrayx);
	Py_CLEAR(arglist);
	//logger("[Callback:R] eval_g");
	return TRUE;
}

Bool
eval_jac_g(Index n, Number * x, Bool new_x,
	   Index m, Index nele_jac,
	   Index * iRow, Index * jCol, Number * values, UserDataPtr data)
{

	//logger("[Callback:E] eval_jac_g");

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	int i;
	long *rowd = NULL;
	long *cold = NULL;

	/* int dims[1]; */
	npy_intp dims[1];
	dims[0] = n;

	double *tempdata;

	if (myowndata->eval_grad_f_python == NULL)	/* Why??? */
		PyErr_Print();

	if (values == NULL) {
		/* import_array (); */
		import_array1(FALSE);

		PyObject *arrayx =
		    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE,
					      (char *)x);
		if (!arrayx)
			return FALSE;

		PyObject *arglist;

		if (user_data != NULL)
			arglist = Py_BuildValue("(OOO)",
						arrayx, Py_True,
						(PyObject *) user_data);
		else
			arglist = Py_BuildValue("(OO)", arrayx, Py_True);

		PyObject *result =
		    PyObject_CallObject(myowndata->eval_jac_g_python, arglist);
		if (!result) {

			logger("[PyIPOPT] return from eval_jac_g is null\n");
			/* TODO: need to deal with reference counting here */
			return FALSE;
		}
		if (!PyTuple_Check(result)) {
			PyErr_Print();
		}
		PyArrayObject *row =
		    (PyArrayObject *) PyTuple_GetItem(result, 0);
		PyArrayObject *col =
		    (PyArrayObject *) PyTuple_GetItem(result, 1);

		if (!row || !col || !PyArray_Check(row) || !PyArray_Check(col)) {
			logger
			    ("[Error] there are problems with row or col in eval_jac_g.\n");
			PyErr_Print();
		}
		rowd = (long *)row->data;
		cold = (long *)col->data;

		for (i = 0; i < nele_jac; i++) {
			iRow[i] = (Index) rowd[i];
			jCol[i] = (Index) cold[i];
		}
		Py_CLEAR(arrayx);
		Py_DECREF(result);
		Py_CLEAR(arglist);
		//logger("[Callback:R] eval_jac_g(1)");
	} else {
		PyObject *arrayx =
		    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE,
					      (char *)x);

		if (!arrayx)
			return FALSE;

		if (new_x && myowndata->apply_new_python) {
			/* Call the python function to applynew */
			PyObject *arg1 = Py_BuildValue("(O)", arrayx);
			PyObject *tempresult =
			    PyObject_CallObject(myowndata->apply_new_python,
						arg1);
			if (tempresult == NULL) {
				logger("[Error] Python function apply_new returns NULL");
				Py_DECREF(arg1);
				return FALSE;
			}
			Py_DECREF(arg1);
			Py_DECREF(tempresult);
		}
		PyObject *arglist;
		if (user_data != NULL)
			arglist = Py_BuildValue("(OOO)",
						arrayx, Py_False,
						(PyObject *) user_data);
		else
			arglist = Py_BuildValue("(OO)", arrayx, Py_False);

		PyArrayObject *result = (PyArrayObject *) PyObject_CallObject(
        myowndata->eval_jac_g_python, arglist);

		if (result == NULL) {
      logger("[Error] Python function eval_jac_g returns NULL");
			PyErr_Print();
      return FALSE;
    }

    if (!PyArray_Check(result)) {
      logger("[Error] Python function eval_jac_g returns non-PyArray");
      Py_DECREF(result);
      return FALSE;
    }

		/*
		 * Code is buggy here. We assume that result is a double
		 * array
		 */
		assert(result->descr->type == 'd');
		tempdata = (double *)result->data;

		for (i = 0; i < nele_jac; i++)
			values[i] = tempdata[i];

		Py_DECREF(result);
		Py_CLEAR(arrayx);
		Py_CLEAR(arglist);
		//logger("[Callback:R] eval_jac_g(2)");
	}
	//logger("[Callback:R] eval_jac_g");
	return TRUE;
}

Bool
eval_h(Index n, Number * x, Bool new_x, Number obj_factor,
       Index m, Number * lambda, Bool new_lambda,
       Index nele_hess, Index * iRow, Index * jCol,
       Number * values, UserDataPtr data)
{
	//logger("[Callback:E] eval_h");

	DispatchData *myowndata = (DispatchData *) data;
	UserDataPtr user_data = (UserDataPtr) myowndata->userdata;

	int i;
	npy_intp dims[1];
	npy_intp dims2[1];

	if (myowndata->eval_h_python == NULL) {
		logger("[Error] There is no eval_h assigned");
		return FALSE;
	}
	if (values == NULL) {
    //logger("[Callback:E] eval_h (1a)");
		PyObject *newx = Py_True;
		PyObject *objfactor = Py_BuildValue("d", obj_factor);
		PyObject *lagrange = Py_True;

		PyObject *arglist;

		if (user_data != NULL) {
			arglist = Py_BuildValue(
          "(OOOOO)", newx, lagrange, objfactor, Py_True,
          (PyObject *) user_data);
    } else {
			arglist = Py_BuildValue(
          "(OOOO)", newx, lagrange, objfactor, Py_True);
    }

    if (arglist == NULL) {
      logger("[Error] failed to build arglist for eval_h");
			PyErr_Print();
      return FALSE;
    } else {
      logger("[Logspam] built arglist for eval_h");
    }

		PyObject *result = PyObject_CallObject(myowndata->eval_h_python, arglist);

    if (result == NULL) {
      logger("[Error] Python function eval_h returns NULL");
			PyErr_Print();
      return FALSE;
    } else {
      logger("[Logspam] Python function eval_h returns non-NULL");
    }

    int result_size = PyTuple_Size(result);

    if (result_size == -1) {
      logger("[Error] Python function eval_h returns non-PyTuple");
      Py_DECREF(result);
      return FALSE;
    }

    if (result_size != 2) {
      logger("[Error] Python function eval_h returns a tuple whose len != 2");
      Py_DECREF(result);
      return FALSE;
    }

    //logger("[Callback:E] eval_h (tuple is the right length)");

		PyArrayObject *row = (PyArrayObject *) PyTuple_GetItem(result, 0);
		PyArrayObject *col = (PyArrayObject *) PyTuple_GetItem(result, 1);

		long *rdata = (long *)row->data;
		long *cdata = (long *)col->data;

		for (i = 0; i < nele_hess; i++) {
			iRow[i] = (Index) rdata[i];
			jCol[i] = (Index) cdata[i];
			/*
			 * logger("PyIPOPT_DEBUG %d, %d\n", iRow[i],
			 * jCol[i]);
			 */
		}

    //logger("[Callback:E] eval_h (clearing stuff now)");

		Py_DECREF(objfactor);
		Py_DECREF(result);
		Py_CLEAR(arglist);
		//logger("[Callback:R] eval_h (1b)");
	} else {
		//logger("[Callback:R] eval_h (2a)");

		PyObject *objfactor = Py_BuildValue("d", obj_factor);

		dims[0] = n;
		PyObject *arrayx =
		    PyArray_SimpleNewFromData(1, dims, PyArray_DOUBLE,
					      (char *)x);
		if (!arrayx)
			return FALSE;

		if (new_x && myowndata->apply_new_python) {
			/* Call the python function to applynew  */
			PyObject *arg1 = Py_BuildValue("(O)", arrayx);
			PyObject *tempresult = PyObject_CallObject(
          myowndata->apply_new_python, arg1);
			if (tempresult == NULL) {
				logger("[Error] Python function apply_new returns NULL");
        PyErr_Print();
				Py_DECREF(arg1);
				return FALSE;
			}
			Py_DECREF(arg1);
			Py_DECREF(tempresult);
		}
		dims2[0] = m;
		PyObject *lagrangex = PyArray_SimpleNewFromData(
        1, dims2, PyArray_DOUBLE, (char *)lambda);
		if (!lagrangex)
			return FALSE;

		PyObject *arglist;

		if (user_data != NULL) {
			arglist = Py_BuildValue(
          "(OOOOO)", arrayx, lagrangex, objfactor, Py_False,
          (PyObject *) user_data);
    } else {
			arglist = Py_BuildValue(
          "(OOOO)", arrayx, lagrangex, objfactor, Py_False);
    }
		PyArrayObject *result = (PyArrayObject *) PyObject_CallObject(
        myowndata->eval_h_python, arglist);

		if (result == NULL) {
      logger("[Error] Python function eval_h returns NULL");
			PyErr_Print();
      return FALSE;
    }

    if (!PyArray_Check(result)) {
      logger("[Error] Python function eval_h returns non-PyArray");
      Py_DECREF(result);
      return FALSE;
    }

		double *tempdata = (double *)result->data;
		for (i = 0; i < nele_hess; i++) {
			values[i] = tempdata[i];
		}
		Py_CLEAR(arrayx);
		Py_CLEAR(lagrangex);
		Py_CLEAR(objfactor);
		Py_DECREF(result);
		Py_CLEAR(arglist);
		//logger("[Callback:R] eval_h (2b)");
	}
	return TRUE;
}
