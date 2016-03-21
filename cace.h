#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#define NS 1

typedef struct {
	int i;
	double v;
} isort;

// cace model object
typedef struct {
    PyObject_HEAD // Python object defaults
	int *col_types;
	PyArrayObject *X;
	PyArrayObject *Y;
	PyArrayObject *ndweights;
	double *weights;
	double *scr[7];
	double *smo;
	double *tx[NS],*ty[NS];
	double *z1,*z2,*z4,*z5;
	double credk;
	isort **xsort;
	int fortran;
	int rows,cols;
	int **card;
} cace;

static void cace_dealloc(cace* self);
static int cace_init(cace *self, PyObject *args, PyObject *kwds);
static PyObject *cace_new(PyTypeObject *type, PyObject *args,PyObject *kwds);
static PyObject *cace_pnew(PyObject *module, PyObject *args);
static PyObject *cace_fit(PyObject *self, PyObject *args,PyObject *keywds);
static PyObject *cace_predict(PyObject *self, PyObject *args,PyObject *keywds);
static PyObject *cace_get_params(PyObject *self, PyObject *args,PyObject *keywds);
static PyObject *cace_set_params(PyObject *self, PyObject *args,PyObject *keywds);

static PyMethodDef cace_methods[] = {
	{"fit",  (PyCFunction)cace_fit, METH_VARARGS|METH_KEYWORDS, "cace fit model to given data."},
	{"predict",  (PyCFunction)cace_predict, METH_VARARGS|METH_KEYWORDS, "cace predict given data."},
	{"get_params",  (PyCFunction)cace_get_params, METH_VARARGS|METH_KEYWORDS, "cace get parameters."},
	{"set_params",  (PyCFunction)cace_set_params, METH_VARARGS|METH_KEYWORDS, "cace set parameters."},
    {NULL}  /* Sentinel */
};

static PyMemberDef cace_members[] = {
    {"credk", T_DOUBLE, offsetof(cace, credk), 0, "cred K"},
/*
    {"weights", T_OBJECT_EX, offsetof(ace, weights), 0, "Training weights"},
    {"offsets", T_OBJECT_EX, offsetof(ace, offsets), 0, "Inital predictions"},
	*/
    {NULL}  /* Sentinel */
};


static PyTypeObject caceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    //0,                         /*ob_size*/
    "cace",             /*tp_name*/
    sizeof(cace),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)cace_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "cace objects",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    cace_methods,             /* tp_methods */
    cace_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)cace_init,      /* tp_init */
    0,                         /* tp_alloc */
    cace_new,                 /* tp_new */
};

