/**
 * @file ace.c - port of Friedman's ACE
 */
#define _CRT_SECURE_NO_WARNINGS
#include "Python.h"
#include <float.h>
#include <math.h>
#include <structmember.h>
#include <omp.h>
#include "cace.h"

#define DELRSQ 0.0001
#define EPS 0.000000001
#define NTERM 3
#define SPAN 0.0
#define ALPHA 5.0
#define CRED_MIN 5.0

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define ISNAN(a) ((!(a<1)) && (!(a>=1)))
#define XC(r,c) (xdata[r*self->cols+c])
#define XVC(r,c) (xvdata[r*self->cols+c])
#define XF(r,c) (xdata[c*self->rows+r])
#define XVF(r,c,rows) (xvdata[c*rows+r])
#define SCC(r,c) (sc[r*5+c])
//#define MC(r,c) (m[c*(self->rows)+r])
#define TXC(is,r,c) (self->tx[is][c*self->rows+r])

double spans[3] = { 0.05, 0.2, 0.5 };

static PyObject *cace_constructor = NULL; 
static void cace_dealloc(cace* self) {
int i;

	Py_XDECREF(self->X);
	Py_XDECREF(self->Y);
	if(self->ndweights) {
		Py_XDECREF(self->ndweights);
	} else {
		if(self->weights) free(self->weights);
	}
	for(i=0;i<7;i++) {
		if(self->scr[i]) {
			free(self->scr[i]);
		}
	}
	if(self->xsort) {
		for(i=0;i<self->cols;i++) free(self->xsort[i]);
		free(self->xsort);
	}
	if(self->card) {
		for(i=0;i<self->cols;i++) 
			if(self->card[i]) free(self->card[i]);
		free(self->card);
	}
	if(self->smo) free(self->smo);
	for(i=0;i<NS;i++) {
		if(self->tx[i]) free(self->tx[i]);
		if(self->ty[i]) free(self->ty[i]);
	}
	if(self->z1) free(self->z1);
	if(self->z2) free(self->z2);
	if(self->z4) free(self->z4);
	if(self->z5) free(self->z5);
	if(self->col_types) free(self->col_types);
    self->ob_base.ob_type->tp_free((PyObject*)self);
}

static PyObject *cace_pickle(PyObject *module, PyObject *args) {

PyObject *obj,*lcol;
int i;
cace *self;

	if (!PyArg_ParseTuple(args, "O!", &caceType, &obj)) return NULL; 
	self=(cace *)obj;
	lcol = Py_BuildValue("[]");
	for(i=0;i<self->cols+1;i++) PyList_Append(lcol,PyLong_FromLong(self->col_types[i]));
	return Py_BuildValue("(O(i,i,i,O,O,O))",cace_constructor,self->rows,self->cols,self->fortran,lcol,self->X,self->Y);
}

//static char *nkwlist[] = {"distribution",NULL };
static PyObject * cace_new(PyTypeObject *type, PyObject *args,PyObject *kwds) {
    cace *self;
	int i;

    self = (cace *)type->tp_alloc(type, 0);
	self->X=NULL;
	self->Y=NULL;
	self->weights=NULL;
	self->ndweights=NULL;
	self->card=NULL;
	self->col_types=NULL;
	self->xsort=NULL;
	for(i=0;i<7;i++) self->scr[i]=NULL;
	for(i=0;i<NS;i++) self->tx[i]=NULL;
	for(i=0;i<NS;i++) self->ty[i]=NULL;
	self->z1=NULL;
	self->z2=NULL;
	self->z4=NULL;
	self->z5=NULL;
	self->smo=NULL;

    return (PyObject *)self;
}
static PyObject *cace_pnew(PyObject *module, PyObject *args) {
cace *self;
PyObject *lcol;
int i,c;

	self = PyObject_NEW(cace, &caceType); 
	if(!PyArg_ParseTuple(args,"|iiiOOO",&self->rows,&self->cols,&self->fortran,&lcol,&self->X,&self->Y)) return NULL;
	self->col_types=(int *)malloc((self->cols+1)*sizeof(int));
	for(c=0;c<self->cols+1;c++) {
		self->col_types[c] = PyLong_AsLong(PyList_GetItem((PyObject *)lcol,c));
	}
	self->X=NULL;
	self->Y=NULL;
	self->weights=NULL;
	self->ndweights=NULL;
	self->col_types=NULL;
	self->xsort=NULL;
	self->card=NULL;
	for(i=0;i<7;i++) self->scr[i]=NULL;
	for(i=0;i<NS;i++) self->tx[i]=NULL;
	for(i=0;i<NS;i++) self->ty[i]=NULL;
	self->z1=NULL;
	self->z2=NULL;
	self->z4=NULL;
	self->z5=NULL;
	self->smo=NULL;
    return (PyObject *)self;
}

static char *ikwlist[] = { NULL };
static int cace_init(cace *self, PyObject *args, PyObject *keywds) {
	if(!PyArg_ParseTupleAndKeywords(args,keywds,"|",ikwlist)) return -1;

    return 0;
}

static PyMethodDef module_methods[] = {
	{"pickle",  (PyCFunction)cace_pickle, METH_VARARGS, "cace pickle support."},
	{"pnew",  (PyCFunction)cace_pnew, METH_VARARGS, "cace unpickle support."},
    {NULL, NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

/*
struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif
*/

#if PY_MAJOR_VERSION >= 3

/*
static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}
*/


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "cace",
        NULL,
        -1, //sizeof(struct module_state),
        module_methods,
        NULL,
        NULL, //myextension_traverse,
        NULL, //myextension_clear,
        NULL
};

#define INITERROR return NULL

PyObject *
PyInit_cace(void)

#else
#define INITERROR return

PyMODINIT_FUNC
initcace(void)
#endif
{
	PyObject *dict;
    PyObject* m;

    if (PyType_Ready(&caceType) < 0) INITERROR;

	#if PY_MAJOR_VERSION >= 3
	    m = PyModule_Create(&moduledef);
	#else
		m = Py_InitModule3("cace", module_methods, "cace fitted model.");
	#endif
    if (m == NULL) INITERROR;

    Py_INCREF(&caceType);
    PyModule_AddObject(m, "cace", (PyObject *)&caceType);
	import_array();

	/*
	PyObject *copy_reg;
	copy_reg = PyImport_ImportModule("copy_reg");
	//cace_constructor=PyObject_GetAttrString(m, "pnew");
	dict = PyModule_GetDict(m); 
	cace_constructor= PyDict_GetItemString(dict, "pnew");
	Py_XINCREF(cace_constructor);
	PyObject_CallMethod(copy_reg,"pickle","OOO",&caceType,PyObject_GetAttrString(m, "pickle"),cace_constructor);
	*/
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}

int cmpi(const void *a,const void *b) {
	isort as = *((isort *)a);
	isort bs = *((isort *)b);
	if(as.v > bs.v) return 1;
	if(as.v < bs.v) return -1;
	return 0;
}

void smooth(cace *self, double *x, double *y, double *w, double span, int iper, double vsmlsq, double *smo, double *avcr) {

double xm, ym, var, cvar, it, tmp, xti, xto, a, h, sy, wt, fbo, fbw;
int i, j, j0, jper, ibw, in, out;

	xm = 0;
	ym = 0;
	var = 0;
	cvar = 0;
	fbw = 0;
	jper = fabs(iper);
	ibw = 0.5 * span * self->rows + 0.5;
	if(ibw < 2) ibw = 2;
	it = 2 * ibw + 1;
	for(i=0; i<it; i++) {
		j = i;
		if(jper == 2) j = i - ibw - 1;
		xti = x[j];
		if(j < 0) {
			j += self->rows;
			xti = x[j] - 1;
		}
		wt = w[j];
		fbo = fbw;
		fbw += wt;
		xm = (fbo * xm + wt * xti) / fbw;
		ym = (fbo * ym + wt * y[j]) / fbw;
		tmp = 0;
		if(fbo > 0.0) tmp = fbw * wt * (xti - xm) / fbo;
		var += tmp * (xti - xm);
		cvar += tmp * (y[j] - ym);
	}
	for(j=0; j<self->rows; j++) {
		out = j - ibw - 1;
		in = j + ibw;
		if(jper==2 || (out >= 0 && in < self->rows)) {
			if(out < 0) {
				out += self->rows;
				xto = x[out] - 1;
				xti = x[in];
			} else {
				if(in >= self->rows) {
					in -= self->rows;
					xti = x[in] + 1;
					xto = x[out];
				} else {
					xti = x[in];
					xto = x[out];
				}
			} 
			wt = w[out];
			fbo = fbw;
			fbw -= wt;
			tmp = 0.0;
			if(fbw > 0) tmp = fbo * wt * (xto - xm) / fbw;
			var -= tmp * (xto - xm);
			cvar -= tmp * (y[out] - ym);
			xm = (fbo * xm - wt * xto) / fbw;
			ym = (fbo * ym - wt * y[out]) / fbw;
			wt = w[in];
			fbo = fbw;
			fbw += wt;
			xm = (fbo * xm + wt * xti) / fbw;
			ym = (fbo * ym + wt * y[in]) / fbw;
			tmp = 0.0;
			if(fbo > 0) tmp = fbw * wt * (xti - xm) / fbo;
			var += tmp * (xti - xm);
			cvar += tmp * (y[in] - ym);
		}
		a = 0;
		if(var > vsmlsq) a = cvar / var;
		smo[j] = a * (x[j] - xm) + ym;
		if(iper > 0) {
			h = 1.0 / fbw;
			if(var > vsmlsq) h += (x[j] - xm)*(x[j] - xm) / var;
			h = min(h,0.9);
			avcr[j] = fabs(y[j] - smo[j]) / (1.0 - w[j] * h);
		}
	}
	j = 0;
	while(j < self->rows) {
		j0 = j;
		sy = smo[j] * w[j];
		fbw = w[j];
		if(j < self->rows-1) {
			while(1) {
				if(x[j+1] > x[j]) break;
				j++;
				sy += w[j] * smo[j];
				fbw += w[j];
				if(j>=self->rows-1) break;
			}
		}
		if(j > j0) {
			sy /= fbw;
			for(i=j0; i<j; i++) {
				smo[i] = sy;
			}
		}
		j++;
	}
}

void supsmu(cace * self, double *x, double *y, double *w, int iper, double span, double *smo) {

double sy,a,scale,*h,vsmlsq,resmin,f,sw;
int i,j,jper;

	if(x[self->rows-1] <= x[0]) {
		sy = 0;
		sw = 0;
		for(j=0;j<self->rows;j++) {
			sy += w[j] * y[j];
			sw += w[j];
		}
		a = sy / sw;
		for(j=0;j<self->rows;j++) {
			smo[j] = a;
		}
		return;
	}

	i = self->rows / 4;
	j = 3 * i;
	scale = x[j] - x[i];
	while(scale <= -1e-200) {
		if(j < self->rows-1) j++;
		if(i >= 0) i--;
		scale = x[j] - x[i];
	}
	vsmlsq = (EPS * scale) * (EPS * scale);
	jper = iper;
	if(iper == 2 && (x[0] < 0.0 || x[self->rows-1] > 1.0)) jper=1;
	if(jper < 1 || jper > 2) jper=1;
	if(span > 0) {
		smooth(self, x, y, w, span, jper, vsmlsq, smo, self->scr[0]);
		return;
	}
	h = (double *)malloc(self->rows * sizeof(double));
	for(i=0; i<3; i++) {
		smooth(self, x, y, w,spans[i], jper, vsmlsq, self->scr[2*i], self->scr[6]);
		smooth(self, x, self->scr[6], w, spans[1], -jper, vsmlsq, self->scr[2*i+1], h);
	}
	for(j=0;j<self->rows;j++) {
		resmin = 1e20;
		for(i=0; i<3; i++) {
			if(self->scr[2*i+1][j] < resmin) {
				resmin = self->scr[2*i+1][j];
				self->scr[6][j] = spans[i];
			}
		}
		if(ALPHA > 0.0 && ALPHA <= 10.0 && resmin < self->scr[5][j]) {
			self->scr[6][j] += (spans[2] - self->scr[6][j]) * pow(max(1e-4, resmin / self->scr[5][j]),10.0-ALPHA);
		}
	}
	smooth(self, x, self->scr[6], w, spans[1], -jper, vsmlsq, self->scr[1], h);
	for(j=0;j<self->rows;j++) {
		if(self->scr[1][j] <= spans[0]) self->scr[1][j] = spans[0];
		if(self->scr[1][j] >= spans[2]) self->scr[1][j] = spans[2];
		f = self->scr[1][j] - spans[1];
		if(f < 0) {
			f = -f / (spans[1] - spans[0]);
			self->scr[3][j] = (1.0 - f) * self->scr[2][j] + f * self->scr[0][j];
		} else {
			f = f / (spans[2] - spans[1]);
			self->scr[3][j] = (1.0 - f) * self->scr[2][j] + f * self->scr[4][j];
		}
	}
	smooth(self, x, self->scr[3], w, spans[0], -jper, vsmlsq, smo, h);
	free(h);
}

void montne(cace *self, double *x) {

int bb, eb, br, er, i, bl, el;
double pmn;

	bb = -1;
	eb = -1;
	while(bb < 1 || x[bb-1] <= x[bb]) {
		if(eb >= self->rows-1) return;
		bb = eb + 1;
		eb = bb;
		while(eb < self->rows-1) {
			if(x[bb] != x[eb+1]) eb++;
		}
LINE30:
		if(eb < self->rows) {
			if(x[eb] > x[eb+1]) {
				br = eb + 1;
				er = br;
				while(er < self->rows-1) {
					if(x[er+1] == x[br]) {
						er++;
					}
				}
				pmn = (x[bb] * (eb - bb + 1) + x[br] * (er - br + 1)) / (er - bb + 1);
				eb = er;
				for(i=bb;i<=eb;i++) {
					x[i] = pmn;
				}
			}
		}
	}
	bl = bb - 1;
	el = bl;
	while(bl >= 1) {
		if(x[bl-1] != x[el]) bl--;
	}
	pmn = (x[bb] * (eb - bb + 1) + x[bl] * (el - bl + 1)) / (eb - bl + 1);
	bb = bl;
	for(i=bb; i<=eb; i++) {
		x[i] = pmn;
	}
	goto LINE30;
}

void smothr(cace *self, int col, int l, double *x, double *y, double *w, double *smo) {

double sm,b,a,sw,d;
int j,j0,i;

	if(l==5) { // category vars
		j = 0;
		while(j < self->rows) {
			j0 = j;
			sm = w[j] * y[j];
			sw = w[j];
			if(j<self->rows-1) {
				while(j<self->rows-1 && x[j+1] <= x[j]) {
					j++;
					sm += w[j]*y[j];
					sw += w[j];
				}
			}
			sm = sm / sw;
			for(i=j0;i<=j;i++) {
				smo[i] = sm;
				if(col>=0) self->card[col][i] = j-j0+1;
			}
			j++;
		}
	} else { // numeric vars
		if(l == 4) { // linear transformation
			sm = 0;
			sw = 0;
			b = 0;
			d = 0;
			for(j=0;j<self->rows;j++) {
				sm += w[j] * x[j] * y[j];
				sw += w[j] * x[j] * x[j];
				b += w[j] * x[j];
				d += w[j];
			}
			a = sm / (sw - (b * b) / d);
			b /= d;
			for(j=0;j<self->rows;j++) {
				smo[j] = a * (x[j] - b);
			}
		} else {
			supsmu(self, x, y, w, l, SPAN, smo);
			if(l == 3) {
				for(j=0;j<self->rows;j++) {
					self->scr[0][j] = smo[j];
					self->scr[1][self->rows-j-1] = smo[j];
				}
				montne(self, self->scr[0]);
				montne(self, self->scr[1]);
				sm = 0;
				sw = 0;
				for(j=0;j<self->rows;j++) {
					sm += (smo[j] - self->scr[0][j]) * (smo[j] - self->scr[0][j]);
					sw += (smo[j] - self->scr[1][self->rows - j - 1]) * (smo[j] - self->scr[1][self->rows - j - 1]);
				}
				if(sm < sw) {
					for(j=0;j<self->rows;j++) {
						smo[j] = self->scr[0][j];
					}
				} else {
					for(j=0;j<self->rows;j++) {
						smo[j] = self->scr[1][self->rows - j - 1];
					}
				}
				j = 0;
				while(j < self->rows) {
					j0 = j;
					sm = smo[j];
					if(j < self->rows) {
						while(x[j+1] <= x[j]) {
							j++;
							sm += smo[j];
							if(j >= self->rows-1) break;
						}
					}
					sm /= j - j0 + 1;
					for(i=j0; i<j; i++) {
						smo[i] = sm;
					}
					j++;
				}
			}
		}
	}
}

void mace(cace *self, int maxit, double *rsq, double *weights) {

int is,i,j,k,iter,js,nt,ism1,start_col;
double *sc,sm,s,sw,*r,ct[NTERM],cmn,cmx,*xdata,*ydata,sv,nit,h,gama,delta,t,u,v,rsqi,sw1;
isort *ysort;
	
	start_col = 0;
RETRY:
	sm = 0;	
	sv = 0;
	sw = 0;
	sw1 = 0;
	xdata=(double *)PyArray_DATA(self->X);
	ydata=(double *)PyArray_DATA(self->Y);
	self->card = (int **)malloc(self->cols*sizeof(int *));
	for(i=0;i<self->cols;i++) {
		self->card[i] = (int *)malloc(self->rows * sizeof(int));
	}
	for(j=0;j<self->rows;j++) {
		sw += weights[j];
	}
	for(is=0;is<NS;is++) {
		for(j=0;j<self->rows;j++) {
			if(self->col_types[self->cols]!=0) {
				self->ty[is][j] = ydata[j];
			}
		}
		for(i=0;i<self->cols;i++) {
			sm=0;
			sw1=0;
			if(self->fortran) {
				for(j=0;j<self->rows;j++) {
					TXC(is,j,i) = XF(j,i);
					sm += weights[j] * XF(j,i);
					sw1 += weights[j];
				}
			} else {
				for(j=0;j<self->rows;j++) {
					TXC(is,j,i) = XC(j,i);
					sm += weights[j] * XC(j,i);
					sw1 += weights[j];
				}
			}
			sm /= sw1;
			for(j=0;j<self->rows;j++) {
				TXC(is,j,i) -= sm;
			}
		}
		sm=0;
		sw1=0;
		for(j=0;j<self->rows;j++) {
			sm += weights[j] * self->ty[is][j];
			sw1 += weights[j];
		}
		sm /= sw1;
		for(j=0;j<self->rows;j++) {
			self->ty[is][j] -= sm;
		}
		sv=0;
		for(j=0;j<self->rows;j++) {
			sv += weights[j]*self->ty[is][j]*self->ty[is][j];
		}
		if(sv==0) {
			PyErr_SetString(PyExc_RuntimeError,"Error: All Y are zero.");
			return;
		}
		if(sv<0.00000001) sv=0.00000001;
		sv = 1/sqrt(sv/sw1);
		for(j=0;j<self->rows;j++) {
			self->ty[is][j] *= sv;
		}

		ysort=(isort *)malloc(sizeof(isort)*self->rows);
		if(is==0) {
			for(j=0;j<self->rows;j++) {
				//MC(j,self->cols) = j;
				//ZC(j,2) = ydata[j];
				ysort[j].i = j;
				ysort[j].v = ydata[j];
			}
			qsort(ysort,self->rows,sizeof(isort),cmpi);
		}
		self->xsort=(isort **)malloc(self->cols*sizeof(isort *));
		for(i=0;i<self->cols;i++) {
			self->xsort[i]=(isort *)malloc(sizeof(isort)*self->rows);
			if(self->fortran) {
				for(j=0;j<self->rows;j++) {
					self->xsort[i][j].i=j;
					self->xsort[i][j].v=XF(j,i);
				}
			} else {
				for(j=0;j<self->rows;j++) {
					self->xsort[i][j].i=j;
					self->xsort[i][j].v=XC(j,i);
				}
			}
			qsort(self->xsort[i],self->rows,sizeof(isort),cmpi);
		}

		// scail
		sc = (double *)malloc(5*self->cols*sizeof(double));
		r = (double *)malloc(self->rows*sizeof(double));
		for(i=0;i<self->cols;i++) {
			SCC(i,0) = 0;
		}
		nit = 0;
		while(1) {
			nit++;
			h=0;
			for(i=0;i<self->cols;i++) {
				SCC(i,4) = SCC(i,0);
				h += SCC(i,0) * SCC(i,0);
			}
			for(iter=0;iter<self->cols;iter++) {
				for(j=0;j<self->rows;j++) {
					s=0;
					for(i=0;i<self->cols;i++) {
						s += SCC(i,0)*TXC(is,j,i);
					}
					r[j]=weights[j]*(self->ty[is][j]-s);
				}
				for(i=0;i<self->cols;i++) {
					s=0;
					for(j=0;j<self->rows;j++) {
						s += r[j]*TXC(is,j,i);
					}
					SCC(i,1) = -2.0*s/sw;
				}
				s = 0;
				for(i=0;i<self->cols;i++) {
					s += SCC(i,1)*SCC(i,1);
				}
				if(s>0) {
					if(iter==0) {
						for(i=0;i<self->cols;i++) {
							SCC(i,2) = -SCC(i,1);
						}
						h = s;
					} else {
						gama = s / h;
						h = s;
						for(i=0;i<self->cols;i++) {
							SCC(i,2) = -SCC(i,1)+gama*SCC(i,3);
						}
					}
					s = 0;
					t = s;
					for(j=0;j<self->rows;j++) {
						u = 0;
						for(i=0;i<self->cols;i++) {
							u = u+SCC(i,2)*TXC(is,j,i);
						}
						s += u*r[j];
						t += weights[j]*u*u;
					}
					delta = s/max(t, 0.00001);
					for(i=0;i<self->cols;i++) {
						SCC(i,0)=SCC(i,0)+delta*SCC(i,2);
						SCC(i,3)=SCC(i,2);
					}
				}
			}
			v = 0;
			for(i=0;i<self->cols;i++) {
				v = max(v,fabs(SCC(i,0)-SCC(i,4)));
			}
			if(v < EPS || nit >= maxit) break;
		} // end while(1)

		for(i=0;i<self->cols;i++) {
			for(j=0;j<self->rows;j++) {
				TXC(is,j,i)=SCC(i,0)*TXC(is,j,i);
			}
		}
		// end scail
		free(sc);
		free(r);

		rsq[is] = 0;
		iter = 0;
		nt = 0;
		for(i=0;i<NTERM;i++) {
			ct[i] = 100.0;
		}
		while(1) {
			iter++;
			nit = 0;
			while(1) {
				rsqi = rsq[is];
				nit++;
				for(j=0;j<self->rows;j++) {
					self->z5[j] = self->ty[is][j];
					for(i=0;i<self->cols;i++) {
						self->z5[j] -= TXC(is,j,i);
					}
				}
				for(i=start_col;i<self->cols;i++) {
					if(self->fortran) {
						for(j=0;j<self->rows;j++) {
							k=self->xsort[i][j].i;
							self->z1[j] = self->z5[k]+TXC(is,k,i);
							self->z2[j] = XF(k,i);
							self->z4[j] = weights[k];
						}
					} else {
						for(j=0;j<self->rows;j++) {
							k=self->xsort[i][j].i;
							self->z1[j] = self->z5[k]+TXC(is,k,i);
							self->z2[j] = XC(k,i);
							self->z4[j] = weights[k];
						}
					}
					smothr(self, i, self->col_types[i], self->z2, self->z1, self->z4, self->smo);
					sm = 0;
					for(j=0;j<self->rows;j++) {
						sm += self->z4[j] * self->smo[j];
					}
					sm /= sw;
					for(j=0;j<self->rows;j++) {
						self->smo[j] -= sm;
					}
					sv = 0;
					for(j=0;j<self->rows;j++) {
						sv += self->z4[j] * (self->z1[j] - self->smo[j]) * (self->z1[j] - self->smo[j]);
					}
					sv = 1 - sv / sw;
					if(sv > rsq[is]) {
						rsq[is] = sv;
						for(j=0;j<self->rows;j++) {
							k = self->xsort[i][j].i;
							TXC(is,k,i) = self->smo[j];
							self->z5[k] = self->z1[j] - self->smo[j];
						}
					}
				}
				if(self->cols == 1 || rsq[is] - rsqi <= DELRSQ || nit >= maxit) break;
			}
			for(j=0;j<self->rows;j++) {
				k=ysort[j].i;
				self->z2[j] = ydata[k];
				self->z4[j] = weights[k];
				self->z1[j] = 0;
				for(i=0;i<self->cols;i++) {
					if(self->col_types[i] != 0) self->z1[j] += TXC(is,k,i);
				}
			}
			smothr(self, -1, self->col_types[self->cols], self->z2, self->z1, self->z4, self->smo);
			if(is > 0) {
				ism1 = is - 1;
				for(js=0; js<ism1; js++) {
					sm = 0;
					for(j=0;j<self->rows;j++) {
						k=ysort[j].i;
						sm += weights[k] * self->smo[j] * self->ty[js][k];
					}
					sm /= sw;
					for(j=0;j<self->rows;j++) {
						k=ysort[j].i;
						self->smo[j] -= sm * self->ty[js][k];
					}
				}
			}
			sm = 0;
			sv = sm;
			for(j=0;j<self->rows;j++) {
				k=ysort[j].i;
				sm += weights[k] * self->smo[j];
				self->z2[k] = self->z1[j];
			}
			sm /= sw;
			for(j=0;j<self->rows;j++) {
				self->smo[j] -= sm;
				sv += self->z4[j] * self->smo[j] * self->smo[j];
			}
			sv /= sw;
			/*
			if(sv <= 0) {
				PyErr_SetString(PyExc_RuntimeError,"Error: sv <= zero.");
				return;
			}
			*/
			if(sv<0.00000001) sv=0.00000001;
			sv = 1 / sqrt(sv);
			for(j=0;j<self->rows;j++) {
				k=ysort[j].i;
				self->ty[is][k] = self->smo[j] * sv;
			}
			sv = 0;
			for(j=0;j<self->rows;j++) {
				sv += weights[j] * (self->ty[is][j] - self->z2[j]) * (self->ty[is][j] - self->z2[j]);
			}
			rsq[is] = 1.0 - (sv / sw);
			nt = nt % NTERM;
			ct[nt] = rsq[is];
			nt++;
			cmn = 100.0;
			cmx = -100.0;
			for(i=0;i<NTERM;i++) {
				cmn = min(cmn, ct[i]);
				cmx = max(cmx, ct[i]);
			}
			if(cmx-cmn <= DELRSQ || iter >= maxit) break;
		} // for iter
		/*
		if(iter == maxit && self->cols == 2) {
			// if 2 columns of data fail to converge, try with just 1
			// this assumes that the first column is less important (ie. an N/A indicator)
			start_col = 1;
			goto RETRY;
		}
		*/
		free(ysort);
	} // for is in NS
}

void model(cace *self, double *f, double *t, double *w) {

int j,i,k,j1,j2;
double s,*ydata;
isort *ysort;

	ydata=(double *)PyArray_DATA(self->Y);
	ysort=(isort *)malloc(sizeof(isort)*self->rows);
	if(self->col_types[self->cols] == 5) {
		for(j=0; j<self->rows; j++) {
			t[j] = self->ty[0][j];
			ysort[j].i = j;
			ysort[j].v = t[j];
		}
	} else {
		for(j=0; j<self->rows; j++) {
			s = 0;
			for(i=0; i<self->cols; i++) {
				s += TXC(0,j,i);
			}
			t[j] = s;
			ysort[j].i = j;
			ysort[j].v = t[j];
		}
	}
	qsort(ysort,self->rows,sizeof(isort),cmpi);
	for(j=0; j<self->rows; j++) {
		t[j] = ysort[j].v;
	}
	for(j=0; j<self->rows; j++) {
		k = ysort[j].i;
		self->z2[j] = w[k];
		if(ydata[k] < 1e20) {
			self->z1[j] = ydata[k];
		} else {
			j1 = j;
			j2 = j;
			while(ydata[ysort[j1].i] > 1e20) {
				j1--;
				if(j1 < 0) break;
			}
			while(ydata[ysort[j2].i] > 1e20) {
				j2++;
				if(j2 >= self->rows) break;
			}
			if(j1 < 0) {
				k = j2;
			} else {
				if(j2 >= self->rows) {
					k = j1;
				} else {
					if(t[j]-t[j1] < t[j2]-t[j]) {
						k = j1;
					} else {
						k = j2;
					}
				}
			}
			self->z1[j] = ydata[ysort[k].i];
			t[j] = t[k];
		}
	}
	if(self->col_types[self->cols] == 5) {
		for(j=0; j<self->rows; j++) {
			f[j] = self->z1[j];
		}
	} else {
		smothr(self, -1, self->col_types[self->cols], t, self->z1, self->z2, f);
	}
	free(ysort);
}

double auto_credk(cace *self, int c) {

double ym,*catvar,*catmean,v,a,*xdata,*ydata;
int maxcat,r,i,*catcount;

	xdata=(double *)PyArray_DATA(self->X);
	ydata=(double *)PyArray_DATA(self->Y);
	ym = 0;
	maxcat=0;
	if(self->fortran) {
		for(r=0;r<self->rows;r++) {
			ym += ydata[r];
			if(XF(r,c)>maxcat) maxcat=XF(r,c);
		}
	} else {
		for(r=0;r<self->rows;r++) {
			ym += ydata[r];
			if(XC(r,c)>maxcat) maxcat=XC(r,c);
		}
	}
	ym /= self->rows;
	catvar = (double *)malloc((maxcat+1)*sizeof(double));
	catmean = (double *)malloc((maxcat+1)*sizeof(double));
	catcount = (int *)malloc((maxcat+1)*sizeof(int));
	for(i=0;i<=maxcat;i++) {
		catmean[i] = 0;
		catvar[i] = 0;
		catcount[i] = 0;
	}
	if(self->fortran) {
		for(r=0;r<self->rows;r++) {
			catmean[(int)XF(r,c)] += ydata[r];
			catcount[(int)XF(r,c)]++;
			catvar[(int)XF(r,c)] += (ydata[r] - ym)*(ydata[r] - ym);
		}
	} else {
		for(r=0;r<self->rows;r++) {
			catmean[(int)XC(r,c)] += ydata[r];
			catcount[(int)XC(r,c)]++;
			catvar[(int)XC(r,c)] += (ydata[r] - ym)*(ydata[r] - ym);
		}
	}
	v=0;
	a=0;
	for(i=0;i<=maxcat;i++) {
		if(catcount[i]) {
			catmean[i] /= catcount[i];
			catvar[i] = catvar[i]/catcount[i];
			v += catcount[i]*catvar[i];
			a += catcount[i]*(catmean[i]-ym)*(catmean[i]-ym);
			//fprintf(stderr,"cat %d %f %f %d %f %f\n",i,catmean[i],catvar[i],catcount[i],v,a);
		}
	}
	v /= self->rows;
	a /= self->rows;
	free(catcount);
	free(catmean);
	free(catvar);
	if(a==0) return CRED_MIN;
	return max(CRED_MIN,v/a);
}

void acemodfc(cace *self, double *f, double *t, int rows, PyArrayObject *X, double *pred, double _credk) {

double th,*xdata,*xvdata,vi,xt,*txm,bc,*credk;
int i,r,sx,place,low,high,jl,jh,rr,card;

	xdata=(double *)PyArray_DATA(self->X);
	xvdata=(double *)PyArray_DATA(X);
	txm = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		txm[i] = 0;
		for(r=0; r<self->rows; r++) {
			txm[i] += TXC(0,r,i);
		}
		txm[i] /= self->rows;
	}
	credk = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		if(_credk == -1 && self->col_types[i]==5)
			credk[i] = auto_credk(self,i);
		else
			credk[i] = _credk;
		self->credk=credk[i];
	}
	for(r=0; r<rows; r++) {
		th = 0;
		for(i=0; i<self->cols; i++) {
			vi = XVC(r,i); // vi = the value we are trying to predict
			if(vi > 1e20) {
				sx = self->xsort[i][self->rows-1].i;
				if(XF(sx,i) > 1e20) th += TXC(0,sx,i);
			} else {
				sx = self->xsort[i][0].i;
				if(vi < XF(sx,i)) {
					place = 0;
					if(self->col_types[i]==5) {
						th += txm[i];
					} else {
						sx = self->xsort[i][place].i;
						th += TXC(0,sx,i);
					}
				} else {
					sx = self->xsort[i][self->rows-1].i;
					if(vi > XF(sx,i)) {
						place = self->rows-1;
						if(self->col_types[i]==5) {
							th += txm[i];
						} else {
							sx = self->xsort[i][place].i;
							th += TXC(0,sx,i);
						}
					} else {
						low = 0;
						high = self->rows-1;
						while(low+1 < high) {
							place = (low + high) / 2;
							sx = self->xsort[i][place].i;
							xt = XF(sx,i);
							if(vi == xt) {
								if(self->col_types[i]==5) {
									card = self->card[i][place];
									bc = card * 1.0 / (card + credk[i]);
									th += bc * TXC(0,sx,i) + (1 - bc) * txm[i];
								} else {
									th += TXC(0,sx,i);
								}
								goto LINE80FC;
							}
							if(vi < xt) {
								high = place;
							} else {
								low = place;
							}
						}
						if(self->col_types[i] != 5) {
							jl = self->xsort[i][low].i;
							jh = self->xsort[i][high].i;
							if(XF(jh,i) > 1e20) {
								th += TXC(0,jl,i);
							} else {
								th += TXC(0,jl,i)+(TXC(0,jh,i)-TXC(0,jl,i))*(vi-XF(jl,i))/(XF(jh,i)-XF(jl,i));
//									  l * (x - l)						     
//							p =	l+h - -----------
//										h - l
							}
						}
					}
				}
			}
LINE80FC:             rr=1;
		} // for i
		if(th <= t[0]) {
			pred[r] = f[0];
			goto NEXTPFC;
		}
		if(th >= t[self->rows-1]) {
			pred[r] = f[self->rows-1];
			goto NEXTPFC;
		}
		low = 0;
		high = self->rows-1;
		while(low+1 < high) {
			place = (low + high) / 2;
			xt = t[place];
			if(th == xt) {
				pred[r] = f[place];
				goto NEXTPFC;
			}
			if(th < xt) {
				high = place;
			} else {
				low = place;
			}
		}
		/*
		if(self->col_types[self->cols] == 5) {
			if(th - t[low] <= t[high] - th) {
				pred[r] = f[low];
				goto NEXTPFC;
			} else {
				pred[r] = f[high];
				goto NEXTPFC;
			}
		}
		*/
		pred[r] = f[low]+(f[high]-f[low])*(th-t[low])/(t[high]-t[low]);
NEXTPFC:  rr=1;
	}
	free(txm);
	free(credk);
}

void acemodcf(cace *self, double *f, double *t, int rows, PyArrayObject *X, double *pred, double _credk) {

double th,*xdata,*xvdata,vi,xt,*txm,bc,*credk;
int i,r,sx,place,low,high,jl,jh,rr,card;

	xdata=(double *)PyArray_DATA(self->X);
	xvdata=(double *)PyArray_DATA(X);
	txm = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		txm[i] = 0;
		for(r=0; r<self->rows; r++) {
			txm[i] += TXC(0,r,i);
		}
		txm[i] /= self->rows;
	}
	credk = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		if(_credk == -1 && self->col_types[i]==5)
			credk[i] = auto_credk(self,i);
		else
			credk[i] = _credk;
		self->credk=credk[i];
	}
	for(r=0; r<rows; r++) {
		th = 0;
		for(i=0; i<self->cols; i++) {
			vi = XVF(r,i,rows); // vi = the value we are trying to predict
			if(vi > 1e20) {
				sx = self->xsort[i][self->rows-1].i;
				if(XC(sx,i) > 1e20) th += TXC(0,sx,i);
			} else {
				sx = self->xsort[i][0].i;
				if(vi < XC(sx,i)) {
					place = 0;
					if(self->col_types[i]==5) {
						th += txm[i];
					} else {
						sx = self->xsort[i][place].i;
						th += TXC(0,sx,i);
					}
				} else {
					sx = self->xsort[i][self->rows-1].i;
					if(vi > XC(sx,i)) {
						place = self->rows-1;
						if(self->col_types[i]==5) {
							th += txm[i];
						} else {
							sx = self->xsort[i][place].i;
							th += TXC(0,sx,i);
						}
					} else {
						low = 0;
						high = self->rows-1;
						while(low+1 < high) {
							place = (low + high) / 2;
							sx = self->xsort[i][place].i;
							xt = XC(sx,i);
							if(vi == xt) {
								if(self->col_types[i]==5) {
									card = self->card[i][place];
									bc = card * 1.0 / (card + credk[i]);
									th += bc * TXC(0,sx,i) + (1 - bc) * txm[i];
								} else {
									th += TXC(0,sx,i);
								}
								goto LINE80CF;
							}
							if(vi < xt) {
								high = place;
							} else {
								low = place;
							}
						}
						if(self->col_types[i] != 5) {
							jl = self->xsort[i][low].i;
							jh = self->xsort[i][high].i;
							if(XC(jh,i) > 1e20) {
								th += TXC(0,jl,i);
							} else {
								th += TXC(0,jl,i)+(TXC(0,jh,i)-TXC(0,jl,i))*(vi-XC(jl,i))/(XC(jh,i)-XC(jl,i));
//									  l * (x - l)						     
//							p =	l+h - -----------
//										h - l
							}
						}
					}
				}
			}
LINE80CF:             rr=1;
		} // for i
		if(th <= t[0]) {
			pred[r] = f[0];
			goto NEXTPCF;
		}
		if(th >= t[self->rows-1]) {
			pred[r] = f[self->rows-1];
			goto NEXTPCF;
		}
		low = 0;
		high = self->rows-1;
		while(low+1 < high) {
			place = (low + high) / 2;
			xt = t[place];
			if(th == xt) {
				pred[r] = f[place];
				goto NEXTPCF;
			}
			if(th < xt) {
				high = place;
			} else {
				low = place;
			}
		}
		/*
		if(self->col_types[self->cols] == 5) {
			if(th - t[low] <= t[high] - th) {
				pred[r] = f[low];
				goto NEXTPCF;
			} else {
				pred[r] = f[high];
				goto NEXTPCF;
			}
		}
		*/
		pred[r] = f[low]+(f[high]-f[low])*(th-t[low])/(t[high]-t[low]);
NEXTPCF:  rr=1;
	}
	free(txm);
	free(credk);
}

void acemodff(cace *self, double *f, double *t, int rows, PyArrayObject *X, double *pred, double _credk) {

double th,*xdata,*xvdata,vi,xt,*txm,bc,*credk;
int i,r,sx,place,low,high,jl,jh,rr,card;

	xdata=(double *)PyArray_DATA(self->X);
	xvdata=(double *)PyArray_DATA(X);
	txm = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		txm[i] = 0;
		for(r=0; r<self->rows; r++) {
			txm[i] += TXC(0,r,i);
		}
		txm[i] /= self->rows;
	}
	credk = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		if(_credk == -1 && self->col_types[i]==5)
			credk[i] = auto_credk(self,i);
		else
			credk[i] = _credk;
		self->credk=credk[i];
	}
	for(r=0; r<rows; r++) {
		th = 0;
		for(i=0; i<self->cols; i++) {
			vi = XVF(r,i,rows); // vi = the value we are trying to predict
			if(vi > 1e20) {
				sx = self->xsort[i][self->rows-1].i;
				if(XF(sx,i) > 1e20) th += TXC(0,sx,i);
			} else {
				sx = self->xsort[i][0].i;
				if(vi < XF(sx,i)) {
					place = 0;
					if(self->col_types[i]==5) {
						th += txm[i];
					} else {
						sx = self->xsort[i][place].i;
						th += TXC(0,sx,i);
					}
				} else {
					sx = self->xsort[i][self->rows-1].i;
					if(vi > XF(sx,i)) {
						place = self->rows-1;
						if(self->col_types[i]==5) {
							th += txm[i];
						} else {
							sx = self->xsort[i][place].i;
							th += TXC(0,sx,i);
						}
					} else {
						low = 0;
						high = self->rows-1;
						while(low+1 < high) {
							place = (low + high) / 2;
							sx = self->xsort[i][place].i;
							xt = XF(sx,i);
							if(vi == xt) {
								if(self->col_types[i]==5) {
									card = self->card[i][place];
									bc = card * 1.0 / (card + credk[i]);
									th += bc * TXC(0,sx,i) + (1 - bc) * txm[i];
								} else {
									th += TXC(0,sx,i);
								}
								goto LINE80FF;
							}
							if(vi < xt) {
								high = place;
							} else {
								low = place;
							}
						}
						if(self->col_types[i] != 5) {
							jl = self->xsort[i][low].i;
							jh = self->xsort[i][high].i;
							if(XF(jh,i) > 1e20) {
								th += TXC(0,jl,i);
							} else {
								th += TXC(0,jl,i)+(TXC(0,jh,i)-TXC(0,jl,i))*(vi-XF(jl,i))/(XF(jh,i)-XF(jl,i));
//									  l * (x - l)						     
//							p =	l+h - -----------
//										h - l
							}
						}
					}
				}
			}
LINE80FF:             rr=1;
		} // for i
		if(th <= t[0]) {
			pred[r] = f[0];
			goto NEXTPFF;
		}
		if(th >= t[self->rows-1]) {
			pred[r] = f[self->rows-1];
			goto NEXTPFF;
		}
		low = 0;
		high = self->rows-1;
		while(low+1 < high) {
			place = (low + high) / 2;
			xt = t[place];
			if(th == xt) {
				pred[r] = f[place];
				goto NEXTPFF;
			}
			if(th < xt) {
				high = place;
			} else {
				low = place;
			}
		}
		/*
		if(self->col_types[self->cols] == 5) {
			if(th - t[low] <= t[high] - th) {
				pred[r] = f[low];
				goto NEXTPFF;
			} else {
				pred[r] = f[high];
				goto NEXTPFF;
			}
		}
		*/
		pred[r] = f[low]+(f[high]-f[low])*(th-t[low])/(t[high]-t[low]);
NEXTPFF:  rr=1;
	}
	free(txm);
	free(credk);
}

void acemodcc(cace *self, double *f, double *t, int rows, PyArrayObject *X, double *pred, double _credk) {

double th,*xdata,*xvdata,vi,xt,*txm,bc,*credk;
int i,r,sx,place,low,high,jl,jh,rr,card;

	xdata=(double *)PyArray_DATA(self->X);
	xvdata=(double *)PyArray_DATA(X);
	txm = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		txm[i] = 0;
		for(r=0; r<self->rows; r++) {
			txm[i] += TXC(0,r,i);
		}
		txm[i] /= self->rows;
	}
	credk = (double *)malloc(self->cols * sizeof(double));
	for(i=0; i<self->cols; i++) {
		if(_credk == -1 && self->col_types[i]==5)
			credk[i] = auto_credk(self,i);
		else
			credk[i] = _credk;
		self->credk=credk[i];
	}
	for(r=0; r<rows; r++) {
		th = 0;
		for(i=0; i<self->cols; i++) {
			vi = XVC(r,i); // vi = the value we are trying to predict
			if(vi > 1e20) {
				sx = self->xsort[i][self->rows-1].i;
				if(XC(sx,i) > 1e20) th += TXC(0,sx,i);
			} else {
				sx = self->xsort[i][0].i;
				if(vi < XC(sx,i)) {
					place = 0;
					if(self->col_types[i]==5) {
						th += txm[i];
					} else {
						sx = self->xsort[i][place].i;
						th += TXC(0,sx,i);
					}
				} else {
					sx = self->xsort[i][self->rows-1].i;
					if(vi > XC(sx,i)) {
						place = self->rows-1;
						if(self->col_types[i]==5) {
							th += txm[i];
						} else {
							sx = self->xsort[i][place].i;
							th += TXC(0,sx,i);
						}
					} else {
						low = 0;
						high = self->rows-1;
						while(low+1 < high) {
							place = (low + high) / 2;
							sx = self->xsort[i][place].i;
							xt = XC(sx,i);
							if(vi == xt) {
								if(self->col_types[i]==5) {
									card = self->card[i][place];
									bc = card * 1.0 / (card + credk[i]);
									th += bc * TXC(0,sx,i) + (1 - bc) * txm[i];
								} else {
									th += TXC(0,sx,i);
								}
								goto LINE80CC;
							}
							if(vi < xt) {
								high = place;
							} else {
								low = place;
							}
						}
						if(self->col_types[i] != 5) {
							jl = self->xsort[i][low].i;
							jh = self->xsort[i][high].i;
							if(XC(jh,i) > 1e20) {
								th += TXC(0,jl,i);
							} else {
								th += TXC(0,jl,i)+(TXC(0,jh,i)-TXC(0,jl,i))*(vi-XC(jl,i))/(XC(jh,i)-XC(jl,i));
//									  l * (x - l)						     
//							p =	l+h - -----------
//										h - l
							}
						}
					}
				}
			}
LINE80CC:             rr=1;
		} // for i
		if(th <= t[0]) {
			pred[r] = f[0];
			goto NEXTPCC;
		}
		if(th >= t[self->rows-1]) {
			pred[r] = f[self->rows-1];
			goto NEXTPCC;
		}
		low = 0;
		high = self->rows-1;
		while(low+1 < high) {
			place = (low + high) / 2;
			xt = t[place];
			if(th == xt) {
				pred[r] = f[place];
				goto NEXTPCC;
			}
			if(th < xt) {
				high = place;
			} else {
				low = place;
			}
		}
		/*
		if(self->col_types[self->cols] == 5) {
			if(th - t[low] <= t[high] - th) {
				pred[r] = f[low];
				goto NEXTPCC;
			} else {
				pred[r] = f[high];
				goto NEXTPCC;
			}
		}
		*/
		pred[r] = f[low]+(f[high]-f[low])*(th-t[low])/(t[high]-t[low]);
NEXTPCC:  rr=1;
	}
	free(txm);
	free(credk);
}


static char *kwlist[] = {"X","y","col_types","maxit","weights",NULL};
static PyObject *cace_fit(PyObject *_self, PyObject *args,PyObject *keywds) {

double *rsq;
cace *self;
char errmsg[256];
int i,c,maxit;
PyListObject *col_types_;
PyObject *rsql;

// Parse Python args into C vars
	self=(cace *)_self;
	Py_XDECREF(self->X);
	Py_XDECREF(self->Y);
	if(self->card) {
		for(i=0;i<self->cols;i++) {
			if(self->card[i]) free(self->card[i]);
		}
		free(self->card);
	}
	if(self->xsort) {
		for(i=0;i<self->cols;i++) free(self->xsort[i]);
		free(self->xsort);
	}
	self->X=NULL;
	self->Y=NULL;
	for(i=0;i<3;i++) {
		if(self->scr[i]) {
			free(self->scr[i]);
			self->scr[i]=NULL;
		}
	}
	for(i=0;i<NS;i++) {
		if(self->tx[i]) {
			free(self->tx[i]);
			self->tx[i]=NULL;
		}
		if(self->ty[i]) {
			free(self->ty[i]);
			self->ty[i]=NULL;
		}
	}
	if(self->smo) {
		free(self->smo);
		self->smo=NULL;
	}
	if(self->z1) {
		free(self->z1);
		self->z1=NULL;
	}
	if(self->z2) {
		free(self->z2);
		self->z2=NULL;
	}
	if(self->z4) {
		free(self->z4);
		self->z4=NULL;
	}
	if(self->z5) {
		free(self->z5);
		self->z5=NULL;
	}
	col_types_=NULL;
	maxit = 20;
	self->ndweights = NULL;
	if(!PyArg_ParseTupleAndKeywords(args,keywds,"O!O!|O!iO",
						kwlist,&PyArray_Type,&self->X,&PyArray_Type,&self->Y,&PyList_Type,&col_types_,
						&maxit,&self->ndweights)) return NULL;
	if(self->X==NULL) {
		PyErr_SetString(PyExc_RuntimeError,"Error: design matrix (X) not specified.");
		return NULL;
	}
	if(self->Y==NULL) {
		PyErr_SetString(PyExc_RuntimeError,"Error: response data (Y) not specified.");
		return NULL;
	}
	if(strcmp(PyArray_DESCR(self->X)->typeobj->tp_name,"numpy.ndarray")) { // X is Numpy ndarray
		PyErr_SetString(PyExc_RuntimeError,"Error: X is not an numpy ndarray");
		return NULL;
	}
	if(strcmp(PyArray_DESCR(self->Y)->typeobj->tp_name,"numpy.ndarray")) { // Y is Numpy ndarray
		PyErr_SetString(PyExc_RuntimeError,"Error: Y is not an numpy ndarray");
		return NULL;
	}
	if(PyArray_NDIM(self->X)!=2) {
		PyErr_SetString(PyExc_RuntimeError,"Error: X does not have two dimensions");
		return NULL;
	}
	if(PyArray_TYPE((PyArrayObject *)self->X)!=NPY_FLOAT64) {
		PyErr_SetString(PyExc_RuntimeError,"Error: X is not an array of 64 bit floats.");
		return NULL;
	}
	if(self->ndweights==Py_None) self->ndweights=NULL;
	if(self->ndweights && strcmp(PyArray_DESCR(self->ndweights)->typeobj->tp_name,"numpy.ndarray")) { // Y is Numpy ndarray
		PyErr_SetString(PyExc_RuntimeError,"Error: weights is not an numpy ndarray");
		return NULL;
	}
	if(self->ndweights && PyArray_TYPE(self->ndweights)!=NPY_FLOAT64) {
		PyErr_SetString(PyExc_RuntimeError,"Error: weights is not an array of 64 bit floats.");
		return NULL;
	}
	self->rows=PyArray_SHAPE(self->X)[0];
	if(self->ndweights && self->rows != PyArray_SHAPE(self->ndweights)[0]) {
		PyErr_SetString(PyExc_RuntimeError,"Error: weights and X have different row counts.");
		return NULL;
	}
	if(self->rows != PyArray_SHAPE(self->Y)[0]) {
		PyErr_SetString(PyExc_RuntimeError,"Error: Y and X have different row counts.");
		return NULL;
	}
	self->cols=PyArray_SHAPE(self->X)[1];
	self->fortran=0;
	if(!(PyArray_FLAGS(self->X) & NPY_ARRAY_C_CONTIGUOUS)) self->fortran=1;
	for(i=0;i<7;i++) {
		self->scr[i]=(double *)malloc(self->rows * sizeof(double));
	}
	for(i=0;i<NS;i++) {
		self->ty[i] = (double *)malloc(self->rows*sizeof(double));
		self->tx[i] = (double *)malloc(self->rows*self->cols*sizeof(double));
	}
	self->z1 = (double *)malloc(self->rows * sizeof(double));
	self->z2 = (double *)malloc(self->rows * sizeof(double));
	self->z4 = (double *)malloc(self->rows * sizeof(double));
	self->z5 = (double *)malloc(self->rows * sizeof(double));
	self->smo=(double *)malloc(self->rows * sizeof(double));

	if(col_types_) {
		if(PyList_Size((PyObject *)col_types_) != self->cols+1) {
			sprintf(errmsg,"Error: length of col_types(%d) != 1 + cols in X(%d)",PyList_Size((PyObject *)col_types_), self->cols+1);
			PyErr_SetString(PyExc_RuntimeError,errmsg);
			return NULL;
		}
	}
	Py_INCREF(self->X);
	Py_INCREF(self->Y);
	if(self->col_types) free(self->col_types);
	self->col_types=(int *)malloc((self->cols+1)*sizeof(int));
	if(col_types_) {
		for(c=0;c<self->cols+1;c++) {
			self->col_types[c] = PyLong_AsLong(PyList_GetItem((PyObject *)col_types_,c));
		}
	} else {
		for(c=0;c<self->cols+1;c++) {
			self->col_types[c] = 1;
		}
	}

	rsq = (double *)malloc(NS*sizeof(double));
	if(self->ndweights) {
		Py_INCREF(self->ndweights);
		self->weights = (double *)PyArray_DATA(self->ndweights);
	} else {
		self->weights = (double *) malloc(self->rows * sizeof(double));
		for(i=0;i<self->rows;i++) {
			self->weights[i]=1;
		}
	}
	mace(self, maxit, rsq, self->weights);
	rsql = Py_BuildValue("[]");
	for(i=0;i<NS;i++) 
		PyList_Append(rsql,Py_BuildValue("d",rsq[i]));
	free(rsq);

	return rsql;
}

static char *tkwlist[] = {"X","credk",NULL};
static PyObject *cace_predict(PyObject *_self, PyObject *args,PyObject *keywds) {

PyArrayObject *X,*Xc;
cace *self;
int rows,cols,fortran;
double *f,*t,*datac,credk;
npy_intp dims[1]; // size of output matrix

// Parse Python args into C vars
	self=(cace *)_self;
	if(self->ty == NULL) {
		PyErr_SetString(PyExc_RuntimeError,"Error: model has not been fit");
		return NULL;
	}
	X=NULL;
	credk = 0;
	if(!PyArg_ParseTupleAndKeywords(args,keywds,"O|d",tkwlist,&X,&credk)) return NULL;
	if(strcmp(PyArray_DESCR(X)->typeobj->tp_name,"numpy.ndarray")) { // X is Numpy ndarray
		PyErr_SetString(PyExc_RuntimeError,"Error: X is not an numpy ndarray");
		return NULL;
	}
	if(PyArray_NDIM(X)!=2) {
		PyErr_SetString(PyExc_RuntimeError,"Error: X does not have two dimensions");
		return NULL;
	}
	if(PyArray_TYPE((PyArrayObject *)X)!=NPY_FLOAT64) {
		PyErr_SetString(PyExc_RuntimeError,"Error: X is not an array of 64 bit floats.");
		return NULL;
	}
	rows=PyArray_SHAPE(X)[0];
	cols=PyArray_SHAPE(X)[1];
	if(cols != self->cols) {
		PyErr_SetString(PyExc_RuntimeError,"Error: predict data has wrong number of cols.");
		return NULL;
	}
	fortran=0;
	if(!(PyArray_FLAGS(X) & NPY_ARRAY_C_CONTIGUOUS)) fortran=1;

	dims[0]=rows;
	Xc=(PyArrayObject *)PyArray_SimpleNew(1,dims,NPY_FLOAT64);
	datac=(double *)PyArray_DATA(Xc);

	f = (double *)malloc(self->rows * sizeof(double));
	t = (double *)malloc(self->rows * sizeof(double));
	model(self, f, t, self->weights);
	/*
	int i;
	double *xdata=(double *)PyArray_DATA(self->X);
	for(i=0;i<self->rows;i++) 
		printf("%f %f %f\n",XC(i,0),f[i],t[i]);
		*/
	if(self->fortran) {
		if(fortran) {
			acemodff(self, f, t, rows, X, datac, credk);
		} else {
			acemodfc(self, f, t, rows, X, datac, credk);
		}
	} else {
		if(fortran) {
			acemodcf(self, f, t, rows, X, datac, credk);
		} else {
			acemodcc(self, f, t, rows, X, datac, credk);
		}
	}
	free(f);
	free(t);

	return (PyObject *)Xc;
}


static char *gkwlist[] = {"deep",NULL };
static PyObject *cace_get_params(PyObject *_self, PyObject *args,PyObject *keywds) {

cace *self;
PyObject *deep;

	self=(cace *)_self;
	if(!PyArg_ParseTupleAndKeywords(args,keywds,"|O!",gkwlist,&PyBool_Type,&deep)) return NULL;
	return Py_BuildValue("");
}
//static char *tkwlist[] = {"distribution",NULL };
static PyObject *cace_set_params(PyObject *_self, PyObject *args,PyObject *keywds) {

cace *self;

	self=(cace *)_self;
	Py_XINCREF(_self);
	return _self;
}
