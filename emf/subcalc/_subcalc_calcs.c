#include <Python.h>
#include <arrayobject.h>
#include "subcalc_calcs.h"

static char module_docstring[] = "Functions for computing the 3-dimensional path integral required to calculate magnetic fields produced by a current carrying wire segment, as part of the Python package emf.subcalc";

static PyObject *subcalc_calcs_compute (PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"compute", subcalc_calcs_compute, METH_VARARGS, NULL}
};

PyMODINIT_FUNC init_compute (void) {
    PyObject *m = Py_InitModule;

    import_array();
}

int main (void) {
    return(1);
}
