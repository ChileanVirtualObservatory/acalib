#include <stdlib.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "morph.h"

static PyObject* morphology_differenceImpl(PyObject* self, PyObject* args);
static PyObject* morphology_segmentationImpl(PyObject* self, PyObject* args);
static PyObject* morphology_erosionImpl(PyObject* self, PyObject* args);

static PyMethodDef module_methods[] =
{
    {"differenceImpl", morphology_differenceImpl, METH_VARARGS, NULL},
    {"segmentationImpl", morphology_segmentationImpl, METH_VARARGS, NULL},
    {"erosionImpl", morphology_erosionImpl, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initmorph()
{
    PyObject* module = Py_InitModule3("morph", module_methods, NULL);
    if(module == NULL)
    {
        return;
    }
    //Load numpy
    import_array();
}

static PyObject* morphology_differenceImpl(PyObject* self, PyObject* args)
{
    PyObject* cumulativeSum_object;
    //Parse input from python
    if(!PyArg_ParseTuple(args, "O", &cumulativeSum_object))
    {
        return NULL;
    }
    //Iterpret input as numpy array
    //https://docs.scipy.org/doc/numpy/reference/c-api.dtype.html#c.NPY_FLOAT64
    PyObject* cumulativeSum_array = PyArray_FROM_OTF(cumulativeSum_object, NPY_FLOAT64, NPY_IN_ARRAY);
    if(cumulativeSum_array == NULL)
    {
        Py_XDECREF(cumulativeSum_array);
        return NULL;
    }
    //Get length of numpy arrray
    int length = (int)PyArray_DIM(cumulativeSum_array, 0);
    //Get pointer to data as C-types
    double* cumulativeSum = (double*)PyArray_DATA(cumulativeSum_array);
    //Call external C function
    double* result = malloc(length*sizeof(double));
    differenceImpl(cumulativeSum, result, length);
    //Clean up
    Py_DECREF(cumulativeSum_array);
    //Create return object
    //https://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html#c.PyArray_SimpleNewFromData
    long dimensions[1];
    dimensions[0] = length;
    PyObject* returnValue = PyArray_SimpleNewFromData(1, dimensions, NPY_FLOAT64, result);
    //Avoid memory leak of the retuned numpy array
    ((PyArrayObject*)returnValue)->flags |= NPY_ARRAY_OWNDATA;
    return returnValue;
}

static PyObject* morphology_segmentationImpl(PyObject* self, PyObject* args)
{
    return Py_BuildValue("");
}

static PyObject* morphology_erosionImpl(PyObject* self, PyObject* args)
{
    return Py_BuildValue("");
}
