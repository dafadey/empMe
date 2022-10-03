#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>

#include <iostream>
#include <vector>
#include <limits>
#include <chrono>

#include "funcall.h"
#include "executor.h"
#include "getGPUtemp.h"

#define USE_SIMPLE_DRAW

#ifdef USE_SIMPLE_DRAW
#include "simpledraw.h"
#endif

#ifndef PYMODNAME
	#define PYMODNAME empMe
#endif

struct empMeObject{
    empMeObject() {
		std::cout << "empMeObject ctor\n";
	}
    PyObject_HEAD
    /* Type-specific fields go here. */
	executor* e_ptr;
	hydro2dHandler* hh_ptr{};
	hydro2dHandler* hh_weak_ptr{};
	int single;
	PyObject* sourceTE{};
	PyObject* sourceTM{};
};

static void
empMe_dealloc(empMeObject *self)
{
  #ifdef USE_SIMPLE_DRAW
  std::cout << "closing simple_draw\n";
  fadey_close();
  #endif
	std::cout << "empMe_dealloc\n";
  if(self->hh_ptr)
    delete self->hh_ptr; // we have to manually delete osc object
  if(self->hh_weak_ptr)
    delete self->hh_weak_ptr; // we have to manually delete osc object
  if(self->e_ptr)
		delete self->e_ptr;
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
empMe_new(PyTypeObject *type, PyObject* args, PyObject* kwds)
{
	std::cout << "empMe_new\n";
  empMeObject *self;
  self = (empMeObject *) type->tp_alloc(type, 0); //NOTE: this does not call ctor just allocates chunk of memory for the structure
  //self = (empMeObject *) type->tp_new(type, args, kwds);
  return (PyObject *) self;
}

double funcTE(double t, double z, void* callerObj) {
	if(!PyCallable_Check(((empMeObject*) callerObj)->sourceTE))
		return PyFloat_AsDouble(0);
	PyObject* res = PyObject_CallFunctionObjArgs(((empMeObject*) callerObj)->sourceTE, PyFloat_FromDouble(t), PyFloat_FromDouble(z), NULL);
	
	double result = PyFloat_AsDouble(res);
	//Py_DECREF(res);
	return result;
}

double funcTM(double t, double z, void* callerObj) {
	//std::cout << "funcTM(" << t << ", " << z << ", " << callerObj << ")\n" << std::flush;
	if(!PyCallable_Check(((empMeObject*) callerObj)->sourceTM))
		return PyFloat_AsDouble(0);
	PyObject* res = PyObject_CallFunctionObjArgs(((empMeObject*) callerObj)->sourceTM, PyFloat_FromDouble(t), PyFloat_FromDouble(z), NULL);

	double result = PyFloat_AsDouble(res);
	//Py_DECREF(res);
	return result;
}

static int
empMe_init(empMeObject *self, PyObject *args, PyObject* kwargs)
{
	std::cout << "empMe_init, self=" << self << "\n";
	static const char* kwlist[] = {"device", "doPhononAbsorbtion", "doPolarization", "extSource", "JHEAT", "single", NULL};
	int device=0;
	int doPhononAbsorbtion = 0;
	int doPolarization = 0;
	int extSource = 0;
	int JHEAT=0;
	self->single = 0;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iiiiii", const_cast<char**>(kwlist), &device, &doPhononAbsorbtion, &doPolarization, &extSource, &JHEAT, &self->single)) {
		std::cerr << "error initalizing hydro handler with given arguments\n";
		return -1;
	}

  self->hh_ptr = new hydro2dHandler(device, static_cast<bool>(doPhononAbsorbtion), static_cast<bool>(doPolarization), static_cast<bool>(extSource), JHEAT);

  self->hh_ptr->TEsource_callback = funcTE;
  self->hh_ptr->TMsource_callback = funcTM;
	self->hh_ptr->pywrapper = (void*) self;
	
	self->e_ptr = new executor(2);
	
  std::cout << "initialized hh with device #" << self->hh_ptr->device << ", PhononAbsorbtion: " << (self->hh_ptr->doPhononAbsorbtion ? "ON" : "OFF") << ", Polarization: " << (self->hh_ptr->doPolarization ? "ON" : "OFF") << ", external Source: " << (self->hh_ptr->extSource ? "ON" : "OFF") << ", JHEAT variant " << self->hh_ptr->JHEAT << ", MODE IS " << (self->single ? "\"SINGLE\"" : "\"NORMAL VS WEAK/LINEAR\"") << '\n';
  
  #ifdef USE_SIMPLE_DRAW
  std::cout << "initializing debugging visual output via simple_draw\n";
  fadey_init(self->hh_ptr->Nx, self->hh_ptr->Nz, 4);
  #endif
  return 0;
}

static PyObject* setup(empMeObject *self, PyObject *args, PyObject* kwargs) {

	std::cout << "empMe_setup, self=" << self << "\n";
	static const char* kwlist[] = {"PMLx","PMLstrength","Lx","surfaceX","Lz","Nx","Nz","switchOnDelay","Bxext","Byext","Bzext","srcTfactor","srcX","srcT","srcTshift","srcNosc","srcNoscTE","srcTTE","srcTshiftTE","srcAperture","srcApertureTE","srcAmp","srcAmpTE","srcPhaseTE","srcPhase","mediaNu","landauDamping","mediaN0","mediaTe0","toothWidth","toothDepth","velocity","mediaDepth","Tavg","NUTratio","cellFilename","cell_scalex","cell_scaley","flip","DRAW","mz_1","mx_1","my_1","diffusion","Tmax","bound_w2","bound_beta","bound_gamma","media_phonon_omega","media_phonon_phw","media_phonon_beta","nonlinVSlin", NULL};


	FL_DBL PMLx;
	FL_DBL Tavg;
	std::string cellFilename;
	FL_DBL cell_scalex;
	FL_DBL cell_scaley;
	int DRAW=0;
	FL_DBL Tmax;
	int nonlinVSlin;
	FL_DBL surfaceX=-1.0;

	#ifdef USEDOUBLE
		const char* argtypes = "|dddddiiddddddddddddddddddddddddddddsddiidddddddddddi";
	#else
		const char* argtypes = "|fffffiiffffffffffffffffffffffffffffsffiifffffffffffi";
	#endif

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, argtypes, const_cast<char**>(kwlist), &PMLx, &self->hh_ptr->PMLstrength, &self->hh_ptr->Lx, &surfaceX, &self->hh_ptr->Lz, &self->hh_ptr->Nx, &self->hh_ptr->Nz, &self->hh_ptr->switchOnDelay, &self->hh_ptr->Bxext, &self->hh_ptr->Byext, &self->hh_ptr->Bzext, &self->hh_ptr->srcTfactor, &self->hh_ptr->srcX, &self->hh_ptr->srcT, &self->hh_ptr->srcTshift, &self->hh_ptr->srcNosc, &self->hh_ptr->srcNoscTE, &self->hh_ptr->srcTTE, &self->hh_ptr->srcTshiftTE, &self->hh_ptr->srcAperture, &self->hh_ptr->srcApertureTE, &self->hh_ptr->srcAmp, &self->hh_ptr->srcAmpTE, &self->hh_ptr->srcPhaseTE, &self->hh_ptr->srcPhase, &self->hh_ptr->mediaNu, &self->hh_ptr->landauDamping, &self->hh_ptr->mediaN0, &self->hh_ptr->mediaTe0, &self->hh_ptr->toothWidth, &self->hh_ptr->toothDepth, &self->hh_ptr->srcVelocity, &self->hh_ptr->mediaDepth, &Tavg, &self->hh_ptr->NUTratio, &cellFilename, &cell_scalex, &cell_scaley, &self->hh_ptr->flip, &DRAW, &self->hh_ptr->mz_1, &self->hh_ptr->mx_1, &self->hh_ptr->my_1, &self->hh_ptr->diffusion, &Tmax, &self->hh_ptr->media_bound_w2, &self->hh_ptr->media_bound_beta, &self->hh_ptr->media_bound_gamma, &self->hh_ptr->media_phonon_omega, &self->hh_ptr->media_phonon_phw, &self->hh_ptr->media_phonon_beta, &nonlinVSlin)) {
		std::cerr << "error initalizing hydro handler with given arguments\n";
		Py_INCREF(Py_None);
			return Py_None;
	}

	self->hh_ptr->dx_1 = FL_DBL(self->hh_ptr->Nx)/self->hh_ptr->Lx;
	self->hh_ptr->dz_1 = FL_DBL(self->hh_ptr->Nz)/self->hh_ptr->Lz;
	self->hh_ptr->dt = FPT(0.5)*((FPT(1.0)/self->hh_ptr->dx_1<FPT(1.0)/self->hh_ptr->dz_1)?FPT(1.0)/self->hh_ptr->dx_1:FPT(1.0)/self->hh_ptr->dz_1);


	if(surfaceX==-1.0)
		surfaceX = self->hh_ptr->Lx / 2.0;

	self->hh_ptr->surfaceX = surfaceX;

	if(self->hh_ptr->srcAmp.size()==0)
		self->hh_ptr->srcAmp.push_back(SRCAMP);
	if(self->hh_ptr->srcAmpTE.size()==0)
		self->hh_ptr->srcAmpTE.push_back(SRCAMP_TE);

	#define FILLTO(REF,ARR, VAR)	for(size_t i=ARR.size(); i < REF.size(); i++)	ARR.push_back(VAR);
	//TM
	FL_DBL lastAperture = self->hh_ptr->srcAperture.size() ? self->hh_ptr->srcAperture[self->hh_ptr->srcAperture.size() - 1] : self->hh_ptr->Lz * 0.25;
	FILLTO(self->hh_ptr->srcAmp, self->hh_ptr->srcAperture, lastAperture);

	FL_DBL lastNosc = self->hh_ptr->srcNosc.size() ? self->hh_ptr->srcNosc[self->hh_ptr->srcNosc.size() - 1] : SRCNOSC;
	FILLTO(self->hh_ptr->srcAmp, self->hh_ptr->srcNosc, lastNosc);

	FILLTO(self->hh_ptr->srcAmp, self->hh_ptr->srcPhase, 0);
	
	FL_DBL lastSrcT = self->hh_ptr->srcT.size() ? self->hh_ptr->srcT[self->hh_ptr->srcT.size() - 1] : SRCT;
	FILLTO(self->hh_ptr->srcAmp, self->hh_ptr->srcT, lastSrcT);

	if(self->hh_ptr->srcTshift.size() != self->hh_ptr->srcT.size()) {
		for(size_t i=self->hh_ptr->srcTshift.size(); i<self->hh_ptr->srcT.size(); i++)
			self->hh_ptr->srcTshift.push_back(self->hh_ptr->srcT[i] * .5);
	}

	//TE
	FL_DBL lastApertureTE = self->hh_ptr->srcApertureTE.size() ? self->hh_ptr->srcApertureTE[self->hh_ptr->srcApertureTE.size() - 1] : self->hh_ptr->Lz * 0.25;
	FILLTO(self->hh_ptr->srcAmp, self->hh_ptr->srcApertureTE, lastApertureTE);
	
	lastNosc = self->hh_ptr->srcNoscTE.size() ? self->hh_ptr->srcNoscTE[self->hh_ptr->srcNoscTE.size() - 1] : SRCNOSC;
	FILLTO(self->hh_ptr->srcAmpTE, self->hh_ptr->srcNoscTE, lastNosc);

	FILLTO(self->hh_ptr->srcAmpTE, self->hh_ptr->srcPhaseTE, 0);

	lastSrcT = self->hh_ptr->srcTTE.size() ? self->hh_ptr->srcTTE[self->hh_ptr->srcTTE.size() - 1] : SRCT;
	FILLTO(self->hh_ptr->srcAmpTE, self->hh_ptr->srcTTE, lastSrcT);

	if(self->hh_ptr->srcTshiftTE.size() != self->hh_ptr->srcT.size()) {
		for(size_t i=self->hh_ptr->srcTshiftTE.size(); i<self->hh_ptr->srcTTE.size(); i++)
			self->hh_ptr->srcTshiftTE.push_back(self->hh_ptr->srcTTE[i]*.5);
	}

	#undef FILLTO


	self->hh_ptr->Bxext=FL_DBL(0);
	self->hh_ptr->Byext=FL_DBL(0);
	self->hh_ptr->Bzext=FL_DBL(0);
	self->hh_ptr->a_cell = nullptr;

	simpleGPUinit(self->hh_ptr);

  if(!self->single)
		self->hh_weak_ptr = new hydro2dHandler(*(self->hh_ptr), self->hh_ptr->device + 1);

	Py_INCREF(Py_None);
  return Py_None;

}

static bool cool_dev(int dev)
{
	int t=getGPUtemp(dev);
	std::cout << "GPU[" << dev << "] temperature=" << t << '\n';
	if(t < 75)
		return true;
	while(getGPUtemp(dev)>75)
		sleep(7);
	printf("temp is %d\n",t);
	return true;
}

static void GPUstep(void* arg) {
	static auto last_time_checked=std::chrono::high_resolution_clock::now();
	auto t=std::chrono::high_resolution_clock::now();
	double dt = (static_cast<std::chrono::duration<double>>(t-last_time_checked)).count();
	if(dt>7) {
		cool_dev(static_cast<hydro2dHandler*>(arg)->device);
		last_time_checked=std::chrono::high_resolution_clock::now();
	}
	simpleGPUstep(static_cast<hydro2dHandler*>(arg));
}

static PyObject* step(empMeObject *self, PyObject** args) {
	static int step_count_global = 0;
	
	self->e_ptr->exec(&GPUstep, (void*) self->hh_ptr, 0);

	if(self->hh_weak_ptr)
		self->e_ptr->exec(&GPUstep, (void*) self->hh_weak_ptr, 1);

	self->e_ptr->sync();

	static std::vector<FL_DBL> ez(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> ex(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> by(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> therm(self->hh_ptr->Nx * self->hh_ptr->Nz);
	if((step_count_global % 100) == 0) {
		std::cout << "step_count_global=" << step_count_global << '\n';
		dev_d2h(self->hh_ptr, self->hh_ptr->Ez, ez.data(), ez.size());
		dev_d2h(self->hh_ptr, self->hh_ptr->Ex, ex.data(), ex.size());
		dev_d2h(self->hh_ptr, self->hh_ptr->By, by.data(), by.size());
		dev_d2h(self->hh_ptr, self->hh_ptr->Te, therm.data(), therm.size());
		
		fadey_draw(ez.data(), self->hh_ptr->Nx, self->hh_ptr->Nz, 0);
		fadey_draw(ex.data(), self->hh_ptr->Nx, self->hh_ptr->Nz, 1);
		fadey_draw(therm.data(), self->hh_ptr->Nx, self->hh_ptr->Nz, 3);
		fadey_draw(by.data(), self->hh_ptr->Nx, self->hh_ptr->Nz, 2);
	}

	step_count_global++;
	
	Py_INCREF(Py_None);
  return Py_None;
}


#define REGSET(X) static PyObject* set_##X(empMeObject *self, PyObject** args) { \
	self->hh_ptr->X = FL_DBL(PyFloat_AsDouble(args[0])); \
	Py_INCREF(Py_None); \
    return Py_None; \
}

#define REGGET(X) static PyObject* get_##X(empMeObject *self) { \
	return PyFloat_FromDouble((double) self->hh_ptr->X); \
}

#define REGGETARR(X) static PyObject* get_array_##X(empMeObject *self, PyObject** args) { \
	int nz = self->hh_ptr->Nz;\
	int nx = self->hh_ptr->Nx;\
	FL_DBL* arr2d = new FL_DBL[nx*nz];\
	dev_d2h(self->hh_ptr, self->hh_ptr->X, arr2d, nx*nz);\
	npy_intp dims[] = {nz,nx};\
	int nd=2;\
	auto newPyArr = PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT32, (void*) arr2d);\
	Py_INCREF(newPyArr);\
	return newPyArr;\
}

#define REGSETARR(X) static PyObject* set_array_##X(empMeObject *self, PyObject** args) { \
}
	//dev_h2d(self->hh_ptr, *args[0], self->hh_ptr->X, self->hh_ptr->Nx * self->hh_ptr->Nz); \

static PyObject* set_surface(empMeObject *self, PyObject* args) {
	std::cout << "hi! set_surface:\n" << std::flush;
	PyObject* pyarr;
	if(PyArg_ParseTuple(args, "O:set_surface", &pyarr)<0)
		std::cout << "parse error\n" << std::flush;
	FL_DBL* arr;
	int dim;
	int n;
	dim = PyArray_NDIM(pyarr);
	std::cout << "array dim is " << dim << '\n';
	std::cout << "PyArray_Type=" << PyArray_TYPE(pyarr) << " NPY_FLOAT32=" << NPY_FLOAT32 << '\n';
	n = PyArray_DIM(pyarr, 0);
	std::cout << "got array from python, size is " << n << '\n' << std::flush;
	PyArray_Descr* descr = PyArray_DESCR(pyarr);
	std::cout << "descr_ptr=" << descr << '\n' << std::flush;
	arr = (FL_DBL*) PyArray_DATA(pyarr);
	FL_DBL min_value = std::numeric_limits<FL_DBL>::max();
	for(int i=0;i<n;i++)
		min_value = std::min(arr[i], min_value);
		
	int nx=self->hh_ptr->Nx;
	int nz=self->hh_ptr->Nz;
	if (nz!= n) {
		std::cout << "ERROR: your array is of wrong dim " << n << ", while Nz=" << nz << '\n';
		Py_INCREF(Py_None);
		return Py_None;
	}
	std::vector<FL_DBL> mat_mask(nx*nz, FL_DBL(0));
	//set zero perturbation for n and Te
	dev_h2d(self->hh_ptr, mat_mask.data(), self->hh_ptr->Te, nx * nz);
	dev_h2d(self->hh_ptr, mat_mask.data(), self->hh_ptr->n, nx * nz);
	
	FL_DBL surfaceX=self->hh_ptr->surfaceX;
	FL_DBL dx=1./self->hh_ptr->dx_1;
	FL_DBL dz=1./self->hh_ptr->dz_1;
	FL_DBL depth=self->hh_ptr->mediaDepth;
	for(int j = self->hh_ptr->PMLimin+3; j<nz-self->hh_ptr->PMLimin-3; j++) {
		for(int i=0;i<nx;i++) {
			FL_DBL x = i * dx;
			if(x>(arr[j]-min_value)+surfaceX && x<(arr[j]-min_value)+surfaceX+depth)
				mat_mask[i+j*nx] = 1.0;
		}
	}
	//set modified material mask
	dev_h2d(self->hh_ptr, mat_mask.data(), self->hh_ptr->mat_mask, nx * nz);

	Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* set_surface_Te(empMeObject *self, PyObject* args) {
	std::cout << "hi! set_surface:\n" << std::flush;
	PyObject* pyarr;
	if(PyArg_ParseTuple(args, "O:set_surface", &pyarr)<0)
		std::cout << "parse error\n" << std::flush;
	FL_DBL* arr;
	int dim;
	int n;
	dim = PyArray_NDIM(pyarr);
	std::cout << "array dim is " << dim << '\n';
	std::cout << "PyArray_Type=" << PyArray_TYPE(pyarr) << " NPY_FLOAT32=" << NPY_FLOAT32 << '\n';
	n = PyArray_DIM(pyarr, 0);
	std::cout << "got array from python, size is " << n << '\n' << std::flush;
	PyArray_Descr* descr = PyArray_DESCR(pyarr);
	std::cout << "descr_ptr=" << descr << '\n' << std::flush;
	arr = (FL_DBL*) PyArray_DATA(pyarr);
	FL_DBL min_value = std::numeric_limits<FL_DBL>::max();
	for(int i=0;i<n;i++)
		min_value = std::min(arr[i], min_value);
		
	int nx=self->hh_ptr->Nx;
	int nz=self->hh_ptr->Nz;
	if (nz!= n) {
		std::cout << "ERROR: your array is of wrong dim " << n << ", while Nz=" << nz << '\n';
		Py_INCREF(Py_None);
		return Py_None;
	}
	std::vector<FL_DBL> mat_mask(nx*nz, FL_DBL(0));
	dev_d2h(self->hh_ptr, self->hh_ptr->mat_mask, mat_mask.data(), nx*nz);
	std::vector<FL_DBL> Te(nx*nz, FL_DBL(0));
	
	FL_DBL surfaceX=self->hh_ptr->surfaceX;
	FL_DBL dx=1./self->hh_ptr->dx_1;
	FL_DBL dz=1./self->hh_ptr->dz_1;
	FL_DBL depth=self->hh_ptr->mediaDepth;
	for(int j = 0; j<nz; j++) {
		for(int i=0;i<nx;i++) {
			FL_DBL x = double(i) * dx;
			if(mat_mask[i+j*nx] == 1.0)
				Te[i+j*nx] = arr[j]*exp(-abs(x-surfaceX)/300.);
		}
	}
	//set modified material mask
	dev_h2d(self->hh_ptr, Te.data(), self->hh_ptr->Te, nx * nz);

	Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* get_therm(empMeObject *self, PyObject** args) {
	int nz = self->hh_ptr->Nz;
	int nx = self->hh_ptr->Nx;
	std::vector<FL_DBL> arr2d(nx*nz, FL_DBL(0));
	FL_DBL* arr = new FL_DBL[nz];
	dev_d2h(self->hh_ptr, self->hh_ptr->Te, arr2d.data(), nx*nz);
	for(int i=0; i<nz; i++) {
		FL_DBL Te=0;
		for(int j=0; j<nx; j++) {
			if(arr2d[nx*i+j]>0) {
				Te += arr2d[nx*i+j];
				break;
			}
		}
		arr[i] = Te;
	}
	npy_intp dims[] = {nz};
	int nd=1;
	auto newPyArr = PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT32, (void*) arr);
	Py_INCREF(newPyArr);
	return newPyArr;
}

REGSET(t)
REGGET(t)
REGSET(dt)
REGGET(dt)
REGGETARR(Ex)
REGSETARR(Ex)
REGGETARR(Te)
REGSETARR(Te)
REGSET(NUTratio)
REGGET(NUTratio)

/*https://docs.python.org/3/c-api/structures.html*/
static PyMemberDef empMeObject_members[] = {
    {"sourceTE", T_OBJECT, offsetof(empMeObject, sourceTE), 0, "source function for TE fields"},
    {"sourceTM", T_OBJECT, offsetof(empMeObject, sourceTM), 0, "source function for TM fields"},
    {NULL}  /* Sentinel */
};

static PyMethodDef empMeObject_methods[] = {
	#define ACCESS(pyname, cname)     {"set_" #pyname, (PyCFunction) set_##cname, METH_FASTCALL, "set data"}, \
  {"get_" #pyname, (PyCFunction) get_##cname, METH_NOARGS, "get data"},
	#define ACCESSARR(pyname, cname)     {"set_2darray_" #pyname, (PyCFunction) set_array_##cname, METH_FASTCALL, "set array"}, \
  {"get_2darray_" #pyname, (PyCFunction) get_array_##cname, METH_NOARGS, "get array"},
    {"step", (PyCFunction) step, METH_FASTCALL, "makes empMe one step"},
    {"setup", (PyCFunction) setup, METH_VARARGS | METH_KEYWORDS, "sets the rig up"},
		ACCESS(t, t)
		ACCESS(dt, dt)
		ACCESS(NUTratio, NUTratio)
		ACCESSARR(Ex, Ex)
		ACCESSARR(Te, Te)
		{"set_surface", (PyCFunction) set_surface, METH_VARARGS, "set surface shape"},
		{"set_surface_Te", (PyCFunction) set_surface_Te, METH_VARARGS, "set surface shape"},
		{"get_therm", (PyCFunction) get_therm, METH_NOARGS, "get integral temperature"},
    {NULL}  /* Sentinel */
    
	#undef ACCESS
	#undef ACCESSARR
};

#define _DEFTOSTR(X) #X
#define DEFTOSTR(X) _DEFTOSTR(X)


static PyTypeObject empMeObjectType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = DEFTOSTR(PYMODNAME) ".hydro2dHandler",
    .tp_basicsize = sizeof(empMeObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor) empMe_dealloc,
		//the following requires initialization of lots of stuff for gcc6. if gcc is updated try to use the below mehtod of initialization and remove unnecessary code from PyInit(void)
    //.tp_flags = Py_TPFLAGS_DEFAULT,
    //.tp_doc = PyDoc_STR("hydrodynamic model with electomagnetic fields"),
    //.tp_methods = empMeObject_methods,
    //.tp_members = empMeObject_members,
    //.tp_init = (initproc) empMe_init,
    //.tp_new = empMe_new,
};

static PyModuleDef empMeModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = DEFTOSTR(PYMODNAME),
    .m_doc = "hydrodynamic model with electomagnetic fields with python controllable external emp source",
    .m_size = -1,
};

#define _MYPYINIT(X) PyMODINIT_FUNC PyInit_ ## X (void)
#define MYPYINIT(X) _MYPYINIT(X)

MYPYINIT(PYMODNAME)
{
		empMeObjectType.tp_flags = Py_TPFLAGS_DEFAULT;
		empMeObjectType.tp_doc = PyDoc_STR("hydrodynamic model with electomagnetic fields");
    empMeObjectType.tp_methods = empMeObject_methods;
    empMeObjectType.tp_members = empMeObject_members;
    empMeObjectType.tp_init = (initproc) empMe_init;
    empMeObjectType.tp_new = empMe_new;
	
    PyObject *m;
    if (PyType_Ready(&empMeObjectType) < 0)
        return NULL;

    m = PyModule_Create(&empMeModule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&empMeObjectType);
    if (PyModule_AddObject(m, "hydro2dHandler", (PyObject *) &empMeObjectType) < 0) {
        Py_DECREF(&empMeObjectType);
        Py_DECREF(m);
        return NULL;
    }
    import_array();

    return m;
}
