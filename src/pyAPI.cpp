#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <chrono>
#include <map>
#include <mutex>

#include "funcall.h"
#include "executor.h"
#include "getGPUtemp.h"
#include "geo.h"

#define USE_SIMPLE_DRAW

#ifdef USE_SIMPLE_DRAW
#include "simpledraw.h"
#endif

#ifndef PYMODNAME
	#define PYMODNAME empMe
#endif


struct timer
{
private:
  std::chrono::high_resolution_clock::time_point t_start;
  std::chrono::high_resolution_clock::time_point t_end;
public:
	double sec() {
		return std::chrono::duration<double>(t_end-t_start).count();
  }
  void start() {
		t_start = std::chrono::high_resolution_clock::now();
  }
  void stop() {
		t_end = std::chrono::high_resolution_clock::now();
  }
};

static PyObject* GPUcount(PyObject* self) {
	return PyLong_FromLong(getGPUcount());
}

static PyObject* GPUtemp(PyObject* self, PyObject** args) {
	int GPUid = PyLong_AsLong(args[0]);
	int temp = getGPUtemp((int)GPUid);
	return PyLong_FromLong((long) temp);
}

static PyObject* GPUtemps(PyObject* self) {
	auto temps = getGPUtemps();
	auto list = PyList_New(temps.size());
	for(int i=0; i < temps.size(); i++)
		PyList_SetItem(list,i,PyLong_FromLong(temps[i]));
	return list;
}

static PyObject* GPUtemps_cool_first(PyObject* self) {
	auto temps = getGPUtemps();
	std::vector<std::array<int,2>> lookup(temps.size());
	for(int i=0; i < temps.size(); i++) {
		lookup[i][0]=temps[i];
		lookup[i][1]=i;
	}
	std::sort(lookup.begin(), lookup.end(), [](const std::array<int,2>& a, const std::array<int,2>& b)->bool {return a[0]<b[0];});
	
	auto list = PyList_New(lookup.size());
	for(int i=0; i < lookup.size(); i++) {
		auto sublist = PyList_New(2);
		PyList_SetItem(sublist,0,PyLong_FromLong(lookup[i][0]));
		PyList_SetItem(sublist,1,PyLong_FromLong(lookup[i][1]));
		PyList_SetItem(list,i,sublist);
	}
	return list;
}

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
	int DRAW;
	PyObject* sourceTE{};
	PyObject* sourceTM{};
	
	std::vector<std::tuple<std::string, FL_DBL**, int>> draw_list;
	std::vector<std::vector<FL_DBL>> host_arrays_collection;	
	
	static int all_draw_lists_size;
	
	int draw_interval;
	int step_count;
	timer tim;
	std::chrono::time_point<std::chrono::high_resolution_clock> last_time_therm_checked;
};

int empMeObject::all_draw_lists_size;

static void
empMe_dealloc(empMeObject *self)
{
  #ifdef USE_SIMPLE_DRAW
  std::cout << "closing simple_draw\n";
  if(self->DRAW)
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
	static std::mutex py_call_mutex;
	py_call_mutex.lock();
	if(!PyCallable_Check(((empMeObject*) callerObj)->sourceTE))
		return PyFloat_AsDouble(0);
	PyObject* res = PyObject_CallFunctionObjArgs(((empMeObject*) callerObj)->sourceTE, PyFloat_FromDouble(t), PyFloat_FromDouble(z), NULL);
	
	double result = PyFloat_AsDouble(res);
	//Py_DECREF(res); // ?
	py_call_mutex.unlock();
	return result;
}

double funcTM(double t, double z, void* callerObj) {
	static std::mutex py_call_mutex;
	py_call_mutex.lock();
	//std::cout << "funcTM(" << t << ", " << z << ", " << callerObj << ")\n" << std::flush;
	if(!PyCallable_Check(((empMeObject*) callerObj)->sourceTM))
		return PyFloat_AsDouble(0);
	PyObject* res = PyObject_CallFunctionObjArgs(((empMeObject*) callerObj)->sourceTM, PyFloat_FromDouble(t), PyFloat_FromDouble(z), NULL);

	double result = PyFloat_AsDouble(res);
	//Py_DECREF(res); // ?
	py_call_mutex.unlock();
	return result;
}

static int
empMe_init(empMeObject *self, PyObject *args, PyObject* kwargs)
{
	std::cout << "empMe_init, self=" << self << "\n";
	static const char* kwlist[] = {"device", "doPhononAbsorbtion", "doPolarization", "extSource", "JHEAT", "single","DRAW", "linear", NULL};
	int device = 0;
	int doPhononAbsorbtion = 0;
	int doPolarization = 0;
	int extSource = 0;
	hydro2dHandler::eHEATTYPE JHEAT=hydro2dHandler::eHEATTYPE::EE;
	self->single = 0;
	self->DRAW = 1;
	self->step_count = 0;
	self->draw_interval = 100;
	empMeObject::all_draw_lists_size = 0;
	
	int islinear=0;
	char* sJHEAT=nullptr;
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iiiisiii", const_cast<char**>(kwlist), &device, &doPhononAbsorbtion, &doPolarization, &extSource, &sJHEAT, &self->single, &self->DRAW, &islinear)) {
		std::cerr << "error initalizing hydro handler with given arguments\n";
		return -1;
	}
  std::cout << "sJHEAT=" << sJHEAT << '\n';
	if(std::string(sJHEAT) == "JJ")
		JHEAT=hydro2dHandler::eHEATTYPE::JJ;
	else if(std::string(sJHEAT) == "JE")
		JHEAT=hydro2dHandler::eHEATTYPE::JE;
	else if(std::string(sJHEAT) == "EE")
		JHEAT=hydro2dHandler::eHEATTYPE::EE;

  self->hh_ptr = new hydro2dHandler(device, static_cast<bool>(doPhononAbsorbtion), static_cast<bool>(doPolarization), static_cast<bool>(extSource), JHEAT);

	self->hh_ptr->TEsource_callback = nullptr;
	self->hh_ptr->TMsource_callback = nullptr;
	self->hh_ptr->pywrapper = (void*) self;
	
	if(islinear)
		self->hh_ptr->linear=true;
	
	self->e_ptr = new executor(2);
	
	self->last_time_therm_checked = std::chrono::high_resolution_clock::now();
	
  std::cout << "initialized hh with device #" << self->hh_ptr->device << ", PhononAbsorbtion: " << (self->hh_ptr->doPhononAbsorbtion ? "ON" : "OFF") << ", Polarization: " << (self->hh_ptr->doPolarization ? "ON" : "OFF") << ", external Source: " << (self->hh_ptr->extSource ? "ON" : "OFF") << ", JHEAT variant " << self->hh_ptr->JHEAT << ", MODE IS " << (self->single ? "\"SINGLE\"" : "\"NORMAL VS WEAK/LINEAR\"") << '\n';
  
  return 0;
}

static PyObject* setup(empMeObject *self, PyObject *args, PyObject* kwargs) {

	std::cout << "empMe_setup, self=" << self << "\n";
	
	static const char* arguments[][2] = {
																			{"PMLx","r"},
																			{"PMLstrength","r"},
																			{"Lx","r"},
																			{"surfaceX","r"},
																			{"Lz","r"},
																			{"Nx","i"},
																			{"Nz","i"},
																			{"switchOnDelay","r"},
																			{"Bxext","r"},
																			{"Byext","r"},
																			{"Bzext","r"},
																			{"srcTfactor","r"},
																			{"srcX","r"},
																			{"srcT","O"},
																			{"srcTshift","O"},
																			{"srcNosc","O"},
																			{"srcNoscTE","O"},
																			{"srcTTE","O"},
																			{"srcTshiftTE","O"},
																			{"srcAperture","O"},
																			{"srcApertureTE","O"},
																			{"srcAmp","O"},
																			{"srcAmpTE","O"},
																			{"srcPhaseTE","O"},
																			{"srcPhase","O"},
																			{"mediaNu","r"},
																			{"landauDamping","r"},
																			{"mediaN0","r"},
																			{"mediaTe0","r"},
																			{"toothWidth","r"},
																			{"toothDepth","r"},
																			{"velocity","r"},
																			{"mediaDepth","r"},
																			{"Tavg","r"},
																			{"NUTratio","r"},
																			{"cellFilename","s"},
																			{"cell_scalex","r"},
																			{"cell_scaley","r"},
																			{"flip","i"},
																			{"mz_1","r"},
																			{"mx_1","r"},
																			{"my_1","r"},
																			{"diffusion","r"},
																			{"Tmax","r"},
																			{"bound_w2","r"},
																			{"bound_beta","r"},
																			{"bound_gamma","r"},
																			{"media_phonon_omega","r"},
																			{"media_phonon_phw","r"},
																			{"media_phonon_beta","r"},
																			{"nonlinVSlin","i"},
																			{"draw_interval","i"}
																			};
	const size_t nargs = sizeof(arguments) / (2 * sizeof(char*));
	static std::vector<const char*> kwlist;
	if(!kwlist.size()) {
	  for(int i=0;i<nargs;i++)
			kwlist.push_back(arguments[i][0]);
		kwlist.push_back(NULL);
	}

	if (self->sourceTE != Py_None)
		self->hh_ptr->TEsource_callback = funcTE;
	if (self->sourceTM != Py_None)
		self->hh_ptr->TMsource_callback = funcTM;

	FL_DBL PMLx;
	FL_DBL Tavg;
	char* cellFilename = nullptr;
	FL_DBL cell_scalex = 1.;
	FL_DBL cell_scaley = 1.;
	FL_DBL Tmax;
	int nonlinVSlin;
	FL_DBL surfaceX=-1.0;

	static std::string argtypes;
	if(argtypes.empty()) {
		argtypes += '|';
	  for(int i=0;i<nargs;i++) {
			const char* v = arguments[i][1];
			if(strcmp(v,"r") == 0)
				#ifdef USEDOUBLE
				argtypes += 'd';
				#else
				argtypes += 'f';
				#endif
			else
				argtypes += v;
		}
	}

	std::cout << "argtypes=\"" << argtypes << "\"\n";
	std::cout << "{";
	for(auto v : kwlist)
			std::cout << "\"" <<  (v ? v : "NULL") << "\", ";
	std::cout << "}\n";

	PyObject* srcT{};
	PyObject* srcTshift{};
	PyObject* srcNosc{};
	PyObject* srcNoscTE{};
	PyObject* srcTTE{};
	PyObject* srcTshiftTE{};
	PyObject* srcAperture{};
	PyObject* srcApertureTE{};
	PyObject* srcAmp{};
	PyObject* srcAmpTE{};
	PyObject* srcPhaseTE{};
	PyObject* srcPhase{};

	if (!PyArg_ParseTupleAndKeywords(args, kwargs, argtypes.c_str(), const_cast<char**>(kwlist.data()), &PMLx, &self->hh_ptr->PMLstrength, &self->hh_ptr->Lx, &surfaceX, &self->hh_ptr->Lz, &self->hh_ptr->Nx, &self->hh_ptr->Nz, &self->hh_ptr->switchOnDelay, &self->hh_ptr->Bxext, &self->hh_ptr->Byext, &self->hh_ptr->Bzext, &self->hh_ptr->srcTfactor, &self->hh_ptr->srcX, &srcT, &srcTshift, &srcNosc, &srcNoscTE, &srcTTE, &srcTshiftTE, &srcAperture, &srcApertureTE, &srcAmp, &srcAmpTE, &srcPhaseTE, &srcPhase, &self->hh_ptr->mediaNu, &self->hh_ptr->landauDamping, &self->hh_ptr->mediaN0, &self->hh_ptr->mediaTe0, &self->hh_ptr->toothWidth, &self->hh_ptr->toothDepth, &self->hh_ptr->srcVelocity, &self->hh_ptr->mediaDepth, &Tavg, &self->hh_ptr->NUTratio, &cellFilename, &cell_scalex, &cell_scaley, &self->hh_ptr->flip, &self->hh_ptr->mz_1, &self->hh_ptr->mx_1, &self->hh_ptr->my_1, &self->hh_ptr->diffusion, &Tmax, &self->hh_ptr->media_bound_w2, &self->hh_ptr->media_bound_beta, &self->hh_ptr->media_bound_gamma, &self->hh_ptr->media_phonon_omega, &self->hh_ptr->media_phonon_phw, &self->hh_ptr->media_phonon_beta, &nonlinVSlin, &self->draw_interval)) {
		std::cerr << "error initalizing hydro handler with given arguments\n" << std::flush;
		Py_INCREF(Py_None);
			return Py_None;
	}

	#define MAKE_INPUT_ARRAY(X) \
															if(X && PyList_Check(X)) { \
																for(int i = 0; i < PyList_Size(X); i++) \
																	self->hh_ptr->X.push_back(FL_DBL(PyFloat_AsDouble(PyList_GetItem(X, i)))); \
															}
	MAKE_INPUT_ARRAY(srcT);
	MAKE_INPUT_ARRAY(srcTshift);
	MAKE_INPUT_ARRAY(srcNosc);
	MAKE_INPUT_ARRAY(srcAperture);
	MAKE_INPUT_ARRAY(srcAmp);
	MAKE_INPUT_ARRAY(srcPhase);
	MAKE_INPUT_ARRAY(srcTTE);
	MAKE_INPUT_ARRAY(srcTshiftTE);
	MAKE_INPUT_ARRAY(srcNoscTE);
	MAKE_INPUT_ARRAY(srcApertureTE);
	MAKE_INPUT_ARRAY(srcAmpTE);
	MAKE_INPUT_ARRAY(srcPhaseTE);

	self->hh_ptr->dx_1 = FL_DBL(self->hh_ptr->Nx)/self->hh_ptr->Lx;
	self->hh_ptr->dz_1 = FL_DBL(self->hh_ptr->Nz)/self->hh_ptr->Lz;
	self->hh_ptr->dt = FPT(0.5)*((FPT(1.0)/self->hh_ptr->dx_1<FPT(1.0)/self->hh_ptr->dz_1)?FPT(1.0)/self->hh_ptr->dx_1:FPT(1.0)/self->hh_ptr->dz_1);

	cell* c_ptr(nullptr);
	if(cellFilename)
	{
		std::cout << "cell file name is " << cellFilename << ", cell_scalex=" << cell_scalex << " ,cell_scaley=" << cell_scaley << ", toothWidth=" << self->hh_ptr->toothWidth << '\n';
		c_ptr = new cell();
		if(!(c_ptr->load_from_svg(cellFilename, self->hh_ptr->toothWidth)))
		{
			std::cout << "Failed to load cel geometry from " << cellFilename << std::endl;
			delete c_ptr;
			c_ptr = nullptr;
		}
		if(c_ptr)
		{
			c_ptr->scalex((double) cell_scalex);
			c_ptr->scaley((double) cell_scaley);
		}
	}
	self->hh_ptr->a_cell=c_ptr;	
	
	if(surfaceX==-1.0)
		surfaceX = self->hh_ptr->Lx / 2.0;

	self->hh_ptr->surfaceX = surfaceX;

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
	simpleGPUstep(static_cast<hydro2dHandler*>(arg));
}

static FL_DBL** field_query(std::string& field_name, hydro2dHandler* hh) {
	static std::map<std::pair<std::string, hydro2dHandler*>, FL_DBL**> lookup;
	auto it = lookup.find(make_pair(field_name, hh));
	if (it == lookup.end()) {
		bool found = false;
		#define REG(x,y) \
		if (field_name == x) { \
			lookup[make_pair(field_name, hh)] = &(hh->y); \
			found = true; \
		}
		
		REG("ex", Ex)
		REG("ey", Ey)
		REG("ez", Ez)
		REG("bx", Bx)
		REG("by", By)
		REG("bz", Bz)
		REG("jx", Jx)
		REG("jy", Jy)
		REG("jz", Jz)
		REG("T", Te)
		REG("n", n)

		#undef REG
		if (found)
			return lookup[make_pair(field_name,hh)];
		else
		  return nullptr;
	} else
		return it->second;
}

static PyObject* sec(empMeObject *self) {
	return PyFloat_FromDouble(self->tim.sec());
}

static PyObject* add_to_drawlist(empMeObject *self, PyObject** args) {
		std::string field_name(PyUnicode_AsUTF8(args[0]));
		FL_DBL** field_pointer = field_query(field_name, self->hh_ptr);
		self->draw_list.push_back(make_tuple(field_name, field_pointer, empMeObject::all_draw_lists_size));
		empMeObject::all_draw_lists_size++;
		Py_INCREF(Py_None);
		return Py_None;
}

static PyObject* step_async(empMeObject *self, PyObject** args) {
	auto t=std::chrono::high_resolution_clock::now();
	double dt = (static_cast<std::chrono::duration<double>>(t - self->last_time_therm_checked)).count();
	if(dt>7) {
		cool_dev(self->hh_ptr->device);
		if(self->hh_weak_ptr)
			cool_dev(self->hh_weak_ptr->device);
		self->last_time_therm_checked=std::chrono::high_resolution_clock::now();
	}

	self->tim.start();
	self->e_ptr->exec(&GPUstep, (void*) self->hh_ptr, 0);

	if(self->hh_weak_ptr)
		self->e_ptr->exec(&GPUstep, (void*) self->hh_weak_ptr, 1);
	
	self->step_count++;
	
	Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* sync_step(empMeObject *self, PyObject** args) {
	
	self->e_ptr->sync();
	self->tim.stop();

  #ifdef USE_SIMPLE_DRAW
  static bool fadey_inited=false;
  if(self->DRAW && !fadey_inited) {
		std::cout << "initializing debugging visual output via simple_draw\n";
		fadey_init(self->hh_ptr->Nx, self->hh_ptr->Nz, empMeObject::all_draw_lists_size);
		fadey_inited=true;
	}
	if((self->step_count % self->draw_interval) == 0 and self->DRAW) {
		if(self->host_arrays_collection.size() != self->draw_list.size())
			self->host_arrays_collection.resize(self->draw_list.size());
		std::cout << "step_count = " << self->step_count << ", drawing: ";
		for(auto draw_item : self->draw_list)
			std::cout << std::get<0>(draw_item) << " ";
		std::cout << '\n';

		for(int i=0; i < self->draw_list.size(); i++) {
			if(!std::get<1>(self->draw_list[i]) || !*std::get<1>(self->draw_list[i]))
				continue;
			if(self->host_arrays_collection[i].size() != self->hh_ptr->Nx * self->hh_ptr->Nz)
				self->host_arrays_collection[i].resize(self->hh_ptr->Nx * self->hh_ptr->Nz);
			dev_d2h(self->hh_ptr, *std::get<1>(self->draw_list[i]), self->host_arrays_collection[i].data(), self->host_arrays_collection[i].size());
			fadey_draw(self->host_arrays_collection[i].data(), self->hh_ptr->Nx, self->hh_ptr->Nz, std::get<2>(self->draw_list[i]));
		}
	}
	
	#endif
	
	Py_INCREF(Py_None);
  return Py_None;
}

static PyObject* step(empMeObject *self, PyObject** args) {

	auto t=std::chrono::high_resolution_clock::now();
	double dt = (static_cast<std::chrono::duration<double>>(t - self->last_time_therm_checked)).count();
	if(dt>7) {
		cool_dev(self->hh_ptr->device);
		if(self->hh_weak_ptr)
			cool_dev(self->hh_weak_ptr->device);
		self->last_time_therm_checked=std::chrono::high_resolution_clock::now();
	}
	
	self->tim.start();
	self->e_ptr->exec(&GPUstep, (void*) self->hh_ptr, 0);

	if(self->hh_weak_ptr)
		self->e_ptr->exec(&GPUstep, (void*) self->hh_weak_ptr, 1);

	self->e_ptr->sync();
	self->tim.stop();

  #ifdef USE_SIMPLE_DRAW
  /*
	static std::vector<FL_DBL> ez(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> ex(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> by(self->hh_ptr->Nx * self->hh_ptr->Nz);
	static std::vector<FL_DBL> therm(self->hh_ptr->Nx * self->hh_ptr->Nz);
	if((self->step_count % 100) == 0 and self->DRAW) {
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
	*/
  static bool fadey_inited=false;
  if(self->DRAW && !fadey_inited) {
		std::cout << "initializing debugging visual output via simple_draw\n";
		fadey_init(self->hh_ptr->Nx, self->hh_ptr->Nz, empMeObject::all_draw_lists_size);
		fadey_inited=true;
	}
	if((self->step_count % self->draw_interval) == 0 and self->DRAW) {
		if(self->host_arrays_collection.size() != self->draw_list.size())
			self->host_arrays_collection.resize(self->draw_list.size());
		std::cout << "step_count = " << self->step_count << ", drawing: ";
		for(auto draw_item : self->draw_list)
			std::cout << std::get<0>(draw_item) << " ";
		std::cout << '\n';

		for(int i=0; i < self->draw_list.size(); i++) {
			if(!std::get<1>(self->draw_list[i]) || !*std::get<1>(self->draw_list[i]))
				continue;
			if(self->host_arrays_collection[i].size() != self->hh_ptr->Nx * self->hh_ptr->Nz)
				self->host_arrays_collection[i].resize(self->hh_ptr->Nx * self->hh_ptr->Nz);
			dev_d2h(self->hh_ptr, *std::get<1>(self->draw_list[i]), self->host_arrays_collection[i].data(), self->host_arrays_collection[i].size());
			fadey_draw(self->host_arrays_collection[i].data(), self->hh_ptr->Nx, self->hh_ptr->Nz, std::get<2>(self->draw_list[i]));
		}
	}
		
	#endif
	self->step_count++;
	
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

#define REGGETi(X) static PyObject* get_##X(empMeObject *self) { \
	return PyLong_FromLong(self->hh_ptr->X); \
}

#define REGGETARR(X) static PyObject* get_array_##X(empMeObject *self, PyObject** args) { \
	int nz = self->hh_ptr->Nz;\
	int nx = self->hh_ptr->Nx;\
	FL_DBL* arr2d = new FL_DBL[nx*nz];\
	dev_d2h(self->hh_ptr, self->hh_ptr->X, arr2d, nx*nz);\
	npy_intp dims[] = {nz,nx};\
	int nd=2;\
	auto newPyArr = PyArray_SimpleNewFromData(nd, dims, NPY_FLOAT32, (void*) arr2d);\
	PyArray_ENABLEFLAGS((PyArrayObject*) newPyArr, NPY_ARRAY_OWNDATA);\
	return newPyArr;\
}
//NOTE:	//Py_INCREF(newPyArr); may be needed above! looks like NO!

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

//fadey
static PyObject* set_surface_Te(empMeObject *self, PyObject* args, PyObject* kwargs) {
	std::cout << "hi! set_surface:\n" << std::flush;
	PyObject* pyarr{};
	FL_DBL depth = 300.;
	const char* kwlist[] = {"surfaceTe", "depth", NULL};
	#ifdef USEDOUBLE
		const char* argtypes = "|Od";
	#else
		const char* argtypes = "|Of";
	#endif
	std::cout << "argtypes=\"" << argtypes << "\"\n";
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, argtypes, const_cast<char**>(kwlist), &pyarr, &depth)) {
		std::cout << "parse error\n" << std::flush;
		return Py_None;
	}
	std::cout << "setting preheat with depth " << depth << '\n';
	FL_DBL* arr;
	int dim;
	int n;
	dim = PyArray_NDIM(pyarr);
	std::cout << "array dim is " << dim << '\n';
	std::cout << "PyArray_Type=" << PyArray_TYPE(pyarr) << ", NPY_FLOAT32=" << NPY_FLOAT32 << ", NPY_FLOAT64=" << NPY_FLOAT64 << '\n';
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
	
	FL_DBL dx=1. / self->hh_ptr->dx_1;
	FL_DBL dz=1. / self->hh_ptr->dz_1;
	for(int j = 0; j<nz; j++) {
		FL_DBL xsurf = .0;
		for(int i=0;i<nx;i++) {
			FL_DBL x = double(i) * dx;
			if(mat_mask[i+j*nx] == 1.0) {
				xsurf = xsurf == .0 ? x : xsurf;
				Te[i+j*nx] = arr[j] * exp( - abs(x - xsurf) / depth);
			}
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
	
	/*
	int nx_surf=0;
	for(; nx_surf<nx && arr2d[nx*int(nz/2)+nx_surf]==0; nx_surf++);

	for(int i=0; i<nz; i++)
		arr[i]=arr2d[nx*i+nx_surf+7];
	*/
	
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
	PyArray_ENABLEFLAGS((PyArrayObject*) newPyArr, NPY_ARRAY_OWNDATA);
	//Py_INCREF(newPyArr); // not needed here, otherwise to much memleak
	return newPyArr;
}

static PyObject* getMaxNu(empMeObject *self) {
	FL_DBL maxNu{0};
	int nz = self->hh_ptr->Nz;
	int nx = self->hh_ptr->Nx;
	static std::vector<FL_DBL> arr2d(nx*nz);
	if(arr2d.size() != nx*nz)
		arr2d.resize(nx*nz);
	dev_d2h(self->hh_ptr, self->hh_ptr->Te, arr2d.data(), nx*nz);
	for(int i=0; i<nx*nz; i++)
		maxNu = std::max(maxNu, self->hh_ptr->mediaNu * (1 + self->hh_ptr->NUTratio * SQRTNUT(arr2d[i])));
	return PyFloat_FromDouble((double) maxNu);
}

static PyObject* get_descr(empMeObject *self) {
	return PyUnicode_FromString(self->hh_ptr->get_description().c_str());
}

#define REG(X) REGSET(X) REGGET(X)
#define REGARR(X) REGSETARR(X) REGGETARR(X)

REGGETi(Nx)
REGGETi(Nz)
REGGET(dx_1)
REGGET(dz_1)
REG(t)
REG(dt)
REG(mediaN0)
REGARR(Jx)
REGARR(Jy)
REGARR(Jz)
REGARR(Ex)
REGARR(Ey)
REGARR(Ez)
REGARR(Bx)
REGARR(By)
REGARR(Bz)
REGARR(Te)
REGARR(n)
REGARR(mat_mask)
REG(NUTratio)
REG(Vdiffusion)
REG(neumann_current_diffusion)

/*https://docs.python.org/3/c-api/structures.html*/
static PyMemberDef empMeObject_members[] = {
    {"sourceTE", T_OBJECT, offsetof(empMeObject, sourceTE), 0, "source function for TE fields"},
    {"sourceTM", T_OBJECT, offsetof(empMeObject, sourceTM), 0, "source function for TM fields"},
    {NULL}  /* Sentinel */
};

static PyMethodDef empMeObject_methods[] = {
	#define ACCESS(pyname, cname)     {"set_" #pyname, (PyCFunction) set_##cname, METH_FASTCALL, "set data"}, \
  {"get_" #pyname, (PyCFunction) get_##cname, METH_NOARGS, "get data"},
	#define ACCESSR(pyname, cname)    {"get_" #pyname, (PyCFunction) get_##cname, METH_NOARGS, "get data"},

	#define ACCESSARR(pyname, cname)     {"set_2darray_" #pyname, (PyCFunction) set_array_##cname, METH_FASTCALL, "set array"}, \
  {"get_2darray_" #pyname, (PyCFunction) get_array_##cname, METH_NOARGS, "get array"},
    {"step", (PyCFunction) step, METH_FASTCALL, "makes empMe one step"},
    {"stepAsync", (PyCFunction) step_async, METH_FASTCALL, "makes empMe one step in asynchronous mode"},
    {"syncStep", (PyCFunction) sync_step, METH_FASTCALL, "performs syncronization and draw"},
    {"setup", (PyCFunction) setup, METH_VARARGS | METH_KEYWORDS, "sets the rig up"},
		{"sec", (PyCFunction) sec, METH_NOARGS, "provides duration of last step in seconds"},
		{"getDescription", (PyCFunction) get_descr, METH_NOARGS, "provides description of all parameters set for hydro2dHandler"},
		ACCESSR(Nx, Nx)
		ACCESSR(Nz, Nz)
		ACCESSR(dx_1, dx_1)
		ACCESSR(dz_1, dz_1)
		ACCESSR(t, t)
		ACCESS(dt, dt)
		ACCESS(NUTratio, NUTratio)
		ACCESS(mediaN0, mediaN0)
		ACCESS(Vdiffusion, Vdiffusion)
		ACCESS(neumann_current_diffusion, neumann_current_diffusion)
		ACCESSARR(Ex, Ex)
		ACCESSARR(Ey, Ey)
		ACCESSARR(Ez, Ez)
		ACCESSARR(Bx, Bx)
		ACCESSARR(By, By)
		ACCESSARR(Bz, Bz)
		ACCESSARR(Jx, Jx)
		ACCESSARR(Jy, Jy)
		ACCESSARR(Jz, Jz)
		ACCESSARR(Te, Te)
		ACCESSARR(n, n)
		ACCESSARR(matMask, mat_mask)
		{"set_surface", (PyCFunction) set_surface, METH_VARARGS, "set surface shape"},
		{"set_surface_Te", (PyCFunction) set_surface_Te, METH_VARARGS | METH_KEYWORDS, "set surface shape"},
		{"get_therm", (PyCFunction) get_therm, METH_NOARGS, "get integral temperature"},
		{"getMaxNu", (PyCFunction) getMaxNu, METH_NOARGS, "get maximum nu to estimate Courants criteria"},
		{"addToDrawList", (PyCFunction) add_to_drawlist, METH_FASTCALL, "adds field to draw list by string name"},
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

static PyMethodDef methods[] = {
	{"GPUcount", (PyCFunction) GPUcount, METH_NOARGS, "provides number of available GPUs"},
	{"GPUtemp", (PyCFunction) GPUtemp, METH_FASTCALL, "provides given GPU temperature"},
	{"GPUtemps", (PyCFunction) GPUtemps, METH_NOARGS, "provides PyList of all the GPU temperatures in system order"},
	{"GPUtempsCoolFirst", (PyCFunction) GPUtemps_cool_first, METH_NOARGS, "provides PyList of all the GPU temperatures with GPU id [[T,id][T,id]...] sorted in accending order with respect to temperature value"},
	{NULL}
};

static PyModuleDef empMeModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = DEFTOSTR(PYMODNAME),
    .m_doc = "hydrodynamic model with electomagnetic fields with python controllable external emp source",
    .m_size = -1,
		methods
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
