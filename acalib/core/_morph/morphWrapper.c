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

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef module_definition =
{
    PyModuleDef_HEAD_INIT,
    "morph",
    NULL,
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_morph()
#else
PyMODINIT_FUNC initmorph()
#endif
{
    #if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&module_definition);
    #else
    PyObject* module = Py_InitModule3("morph", module_methods, NULL);
    #endif
    if(module == NULL)
    {
        return NULL;
    }
    import_array();
    #if PY_MAJOR_VERSION >= 3
    return module;
    #endif
}

static PyObject* morphology_differenceImpl(PyObject* self, PyObject* args)
{
    PyObject* cumulativeSum_object;
    if(!PyArg_ParseTuple(args, "O", &cumulativeSum_object))
    {
        return NULL;
    }
    PyObject* cumulativeSum_array = PyArray_FROM_OTF(cumulativeSum_object, NPY_FLOAT64, NPY_IN_ARRAY);
    if(cumulativeSum_array == NULL)
    {
        Py_XDECREF(cumulativeSum_array);
        return NULL;
    }
    int length = (int)PyArray_DIM(cumulativeSum_array, 0);
    double* cumulativeSum = (double*)PyArray_DATA(cumulativeSum_array);
    double* result = malloc(length*sizeof(double));
    differenceImpl(cumulativeSum, result, length);
    Py_DECREF(cumulativeSum_array);
    long dimensions[1];
    dimensions[0] = length;
    PyObject* returnValue = PyArray_SimpleNewFromData(1, dimensions, NPY_FLOAT64, result);
    ((PyArrayObject*)returnValue)->flags |= NPY_ARRAY_OWNDATA;
    return returnValue;
}

static PyObject* morphology_segmentationImpl(PyObject* self, PyObject* args)
{
    PyObject* input;
    if(!PyArg_ParseTuple(args, "O", &input))
    {
        return NULL;
    }
    PyObject* input_array = PyArray_FROM_OTF(input, NPY_FLOAT64, NPY_IN_ARRAY);
    if(input_array == NULL)
    {
        Py_XDECREF(input_array);
        return NULL;
    }
    int length = (int)PyArray_DIM(input_array, 0);
    double* input_data = (double*)PyArray_DATA(input_array);
    double* result = malloc(length*sizeof(double));
    segmentationImpl(input_data, result, length);
    Py_DECREF(input_array);
    long dimensions[1];
    dimensions[0] = length;
    PyObject* returnValue = PyArray_SimpleNewFromData(1, dimensions, NPY_FLOAT64, result);
    ((PyArrayObject*)returnValue)->flags |= NPY_ARRAY_OWNDATA;
    return returnValue;
}

static PyObject* morphology_erosionImpl(PyObject* self, PyObject* args)
{
    PyObject* input;
    if(!PyArg_ParseTuple(args, "O", &input))
    {
        return NULL;
    }
    PyObject* input_array = PyArray_FROM_OTF(input, NPY_FLOAT64, NPY_IN_ARRAY);
    if(input_array == NULL)
    {
        Py_XDECREF(input_array);
        return NULL;
    }
    int length = (int)PyArray_DIM(input_array, 0);
    double* input_data = (double*)PyArray_DATA(input_array);
    double* result = malloc(length*sizeof(double));
    erosionImpl(input_data, result, length);
    memcpy(result, input_data, length*sizeof(double));
    Py_DECREF(input_array);
    long dimensions[1];
    dimensions[0] = length;
    PyObject* returnValue = PyArray_SimpleNewFromData(1, dimensions, NPY_FLOAT64, result);
    ((PyArrayObject*)returnValue)->flags |= NPY_ARRAY_OWNDATA;
    return returnValue;
}
