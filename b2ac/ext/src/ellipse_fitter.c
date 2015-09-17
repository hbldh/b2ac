/**
 **************************************************************************************
 * @file    $ellipse_fitter.c
 * @author  Henrik Blidh
 * @version 1.0
 * @date    2013-08-01
 * @brief   Numpy/C extensions for using ellipse fitting methods.
 **************************************************************************************
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL nptest_ARRAY_API
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "EllipseFit.h"
#include "EllipseFit_double.h"
#include "EllipseFit_float.h"
#include "EllipseFit_int.h"

static char fit_ellipse_double_docs[] =
      "Algorithm for fitting ellipses.\n\n"
      "Definition:\n"
      "  ``fit_ellipse_double(points)``\n\n"
      ":param points: X- and y-coordinate array, size [N x 2], of the points to fit ellipse to.\n"
      ":type points: :py:class:`numpy.ndarray`\n"
      ":returns: Center point, x and y-axis semiaxes and the rotation angle.\n"
      ":rtype: tuple\n\n";

static PyObject *fit_ellipse_double_func(PyObject *self, PyObject *args)
{
    PyObject *points=NULL;
    PyObject *points_arr=NULL;

    /* Local variables for handling the dimensions of the ndarrays. */
    long *points_arr_dims=NULL;
    int points_arr_nd=0;

    /* Pointers into the ndarrays. */
    npy_int32 *points_arr_ptr=NULL;

    /* Loop counters and other local variables. */
    int i;

    /* Parse the input to correct pointers. */
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;

    points_arr = PyArray_FROM_OTF(points, NPY_INT32, NPY_IN_ARRAY);
    if (points_arr == NULL) goto fail;

    /* Dimensions of the arrays. */
    points_arr_nd        = PyArray_NDIM(points_arr);
    points_arr_dims      = PyArray_DIMS(points_arr);

    if (points_arr_nd != 2)
    {
        printf("Dimension of coordinates array must be equal to zero. (%d != 2).\n",
               points_arr_nd);
        goto fail;
    }

    int32_t *x = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    int32_t *y = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    for (i = 0; i < points_arr_dims[0]; i++)
    {
        points_arr_ptr = PyArray_GETPTR2(points_arr, i, 0);
        x[i] = points_arr_ptr[0];
        y[i] = points_arr_ptr[1];
    }

    // Call the actual ellipse fitting function.
    EllipseDouble ellipse = FitEllipse_double(x, y, (int16_t) points_arr_dims[0]);

    /* Decrease reference counters and free memory to avoid memory leaks. */
    free(x);
    free(y);
    Py_DECREF(points_arr);

    if (ellipse.isValid)
    {
        int out_dims[1] = {2};
        PyArrayObject *xy;
        npy_float64 *xy_ptr;
        xy = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        xy_ptr = (npy_float64 *)PyArray_GETPTR1(xy, 0);
        xy_ptr[0] = ellipse.centerX;
        xy_ptr[1] = ellipse.centerY;

        PyArrayObject *axes;
        npy_float64 *axes_ptr;
        axes = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        axes_ptr = (npy_float64 *)PyArray_GETPTR1(axes, 0);
        axes_ptr[0] = ellipse.xAxis;
        axes_ptr[1] = ellipse.yAxis;

        PyObject *tuple_result = PyTuple_New(3);
        PyTuple_SetItem(tuple_result, 0, PyArray_Return(xy));
        PyTuple_SetItem(tuple_result, 1, PyArray_Return(axes));
        PyTuple_SetItem(tuple_result, 2, PyFloat_FromDouble(ellipse.rotationAngle));

        return (PyObject *)tuple_result;
    }
    else
    {
        Py_RETURN_NONE;
    }

 fail:
    /* If error occurs, decrease reference counters to avoid memory leaks. */
    Py_XDECREF(points_arr);
    /* Return error indication. */
    return NULL;
}

static char fit_ellipse_float_docs[] =
      "Algorithm for fitting ellipses.\n\n"
      "Definition:\n"
      "  ``fit_ellipse_float(points)``\n\n"
      ":param points: X- and y-coordinate array, size [N x 2], of the points to fit ellipse to.\n"
      ":type points: :py:class:`numpy.ndarray`\n"
      ":returns: Center point, x and y-axis semiaxes and the rotation angle.\n"
      ":rtype: tuple\n\n";

static PyObject *fit_ellipse_float_func(PyObject *self, PyObject *args)
{
    PyObject *points=NULL;
    PyObject *points_arr=NULL;

    /* Local variables for handling the dimensions of the ndarrays. */
    long *points_arr_dims=NULL;
    int points_arr_nd=0;

    /* Pointers into the ndarrays. */
    npy_int32 *points_arr_ptr=NULL;

    /* Loop counters and other local variables. */
    int i;

    /* Parse the input to correct pointers. */
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;

    points_arr = PyArray_FROM_OTF(points, NPY_INT32, NPY_IN_ARRAY);
    if (points_arr == NULL) goto fail;

    /* Dimensions of the arrays. */
    points_arr_nd        = PyArray_NDIM(points_arr);
    points_arr_dims      = PyArray_DIMS(points_arr);

    if (points_arr_nd != 2)
    {
        printf("Dimension of coordinates array must be equal to zero. (%d != 2).\n",
               points_arr_nd);
        goto fail;
    }

    int32_t *x = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    int32_t *y = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    for (i = 0; i < points_arr_dims[0]; i++)
    {
        points_arr_ptr = PyArray_GETPTR2(points_arr, i, 0);
        x[i] = points_arr_ptr[0];
        y[i] = points_arr_ptr[1];
    }

    // Call the actual ellipse fitting function.
    EllipseFloat ellipse = FitEllipse_float(x, y, (int16_t) points_arr_dims[0]);

    /* Decrease reference counters and free memory to avoid memory leaks. */
    free(x);
    free(y);
    Py_DECREF(points_arr);

    if (ellipse.isValid)
    {
        int out_dims[1] = {2};
        PyArrayObject *xy;
        npy_float64 *xy_ptr;
        xy = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        xy_ptr = (npy_float64 *)PyArray_GETPTR1(xy, 0);
        xy_ptr[0] = ellipse.centerX;
        xy_ptr[1] = ellipse.centerY;

        PyArrayObject *axes;
        npy_float64 *axes_ptr;
        axes = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        axes_ptr = (npy_float64 *)PyArray_GETPTR1(axes, 0);
        axes_ptr[0] = ellipse.xAxis;
        axes_ptr[1] = ellipse.yAxis;

        PyObject *tuple_result = PyTuple_New(3);
        PyTuple_SetItem(tuple_result, 0, PyArray_Return(xy));
        PyTuple_SetItem(tuple_result, 1, PyArray_Return(axes));
        PyTuple_SetItem(tuple_result, 2, PyFloat_FromDouble(ellipse.rotationAngle));

        return (PyObject *)tuple_result;
    }
    else
    {
        Py_RETURN_NONE;
    }

 fail:
    /* If error occurs, decrease reference counters to avoid memory leaks. */
    Py_XDECREF(points_arr);
    /* Return error indication. */
    return NULL;
}

static char fit_ellipse_int_docs[] =
      "Algorithm for fitting ellipses.\n\n"
      "Definition:\n"
      "  ``fit_ellipse_int(points)``\n\n"
      ":param points: X- and y-coordinate array, size [N x 2], of the points to fit ellipse to.\n"
      ":type points: :py:class:`numpy.ndarray`\n"
      ":returns: Center point, x and y-axis semiaxes and the rotation angle.\n"
      ":rtype: tuple\n\n";

static PyObject *fit_ellipse_int_func(PyObject *self, PyObject *args)
{
    PyObject *points=NULL;
    PyObject *points_arr=NULL;

    /* Local variables for handling the dimensions of the ndarrays. */
    long *points_arr_dims=NULL;
    int points_arr_nd=0;

    /* Pointers into the ndarrays. */
    npy_int32 *points_arr_ptr=NULL;

    /* Loop counters and other local variables. */
    int i;

    /* Parse the input to correct pointers. */
    if (!PyArg_ParseTuple(args, "O", &points))
        return NULL;

    points_arr = PyArray_FROM_OTF(points, NPY_INT32, NPY_IN_ARRAY);
    if (points_arr == NULL) goto fail;

    /* Dimensions of the arrays. */
    points_arr_nd        = PyArray_NDIM(points_arr);
    points_arr_dims      = PyArray_DIMS(points_arr);

    if (points_arr_nd != 2)
    {
        printf("Dimension of coordinates array must be equal to zero. (%d != 2).\n",
               points_arr_nd);
        goto fail;
    }

    int32_t *x = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    int32_t *y = (int32_t*) malloc(sizeof(int32_t) * points_arr_dims[0]);
    for (i = 0; i < points_arr_dims[0]; i++)
    {
        points_arr_ptr = PyArray_GETPTR2(points_arr, i, 0);
        x[i] = points_arr_ptr[0];
        y[i] = points_arr_ptr[1];
    }

    // Call the actual ellipse fitting function.
    EllipseFloat ellipse = FitEllipse_int(x, y, (int16_t) points_arr_dims[0]);

    /* Decrease reference counters and free memory to avoid memory leaks. */
    free(x);
    free(y);
    Py_DECREF(points_arr);

    if (ellipse.isValid)
    {
        int out_dims[1] = {2};
        PyArrayObject *xy;
        npy_float64 *xy_ptr;
        xy = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        xy_ptr = (npy_float64 *)PyArray_GETPTR1(xy, 0);
        xy_ptr[0] = ellipse.centerX;
        xy_ptr[1] = ellipse.centerY;

        PyArrayObject *axes;
        npy_float64 *axes_ptr;
        axes = (PyArrayObject *) PyArray_FromDims(1, out_dims, NPY_FLOAT64);
        axes_ptr = (npy_float64 *)PyArray_GETPTR1(axes, 0);
        axes_ptr[0] = ellipse.xAxis;
        axes_ptr[1] = ellipse.yAxis;

        PyObject *tuple_result = PyTuple_New(3);
        PyTuple_SetItem(tuple_result, 0, PyArray_Return(xy));
        PyTuple_SetItem(tuple_result, 1, PyArray_Return(axes));
        PyTuple_SetItem(tuple_result, 2, PyFloat_FromDouble(ellipse.rotationAngle));

        return (PyObject *)tuple_result;
    }
    else
    {
        Py_RETURN_NONE;
    }

 fail:
    /* If error occurs, decrease reference counters to avoid memory leaks. */
    Py_XDECREF(points_arr);
    /* Return error indication. */
    return NULL;
}


static PyMethodDef EFMethods[] = {
    {"fit_ellipse_double", fit_ellipse_double_func, METH_VARARGS, fit_ellipse_double_docs},
    {"fit_ellipse_float", fit_ellipse_float_func, METH_VARARGS, fit_ellipse_float_docs},
    {"fit_ellipse_int", fit_ellipse_int_func, METH_VARARGS, fit_ellipse_int_docs},
     {NULL, NULL, 0, NULL} /* Sentinel */
};

PyMODINIT_FUNC initellipse_fitter(void)
{
    char * docstring = "Ellipse Fitting extension.";
    (void) Py_InitModule3("ellipse_fitter", EFMethods, docstring);
    import_array();
}
