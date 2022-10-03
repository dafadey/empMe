TARGET_bin := test_mGPU
TARGET_STATIC_bin := test_mGPU_static
TARGET_lib := pyEmpMe.so
TARGET_STATIC_lib := pyEmpMe_static.so

arch :=CC30 #CC35
mrch :=sm_30 #sm_35
enable_double := 
MYTOOLS := /home/dan/tools
INC := -I/usr/local/cuda-9.0/include -I/usr/local/cuda-9.0/samples/common/inc -I$(MYTOOLS)/include -I$(MYTOOLS)/include/mesh1 #-I$(MYTOOLS)/nvml/include 
cuLIB := /usr/local/cuda-9.0/lib64
LIB := -L$(cuLIB) -L$(MYTOOLS)/lib #-L$(MYTOOLS)/nvml/lib64
fm := #--use_fast_math #coment it if you whant divisions and transcend funcs to be IEEE or something similar to IEEE
MYLIBS := -lsimpledraw_glfw3 -lzstream -lsvg1 -lmesh1 -ltclient -lgetGPUtemp

gccOPTS := -fPIC -std=c++11 -O3

NVCC_COMMAND_LINE := /usr/local/cuda-9.0/bin/nvcc -w -Xcompiler -fPIC $(enable_double) $(INC) -O3 -c -arch=$(mrch) $(fm)

GCC_COMMAND_LINE := /usr/bin/g++-6 -g $(INC) $(gccOPTS) $(enable_double) -Wall -Werror -c

GCC_PY_COMMAND_LINE := /usr/bin/g++-6 -g $(INC) $(gccOPTS) -I/usr/include/python3.9 -c

LINK_COMMAND_LINE := /usr/bin/g++-6 -fPIC $(LIB) $(MYLIBS) -Wl,-rpath=$(MYTOOLS)/lib -Wl,-rpath=${cuLIB} -fopenmp -lpthread -lcudart

.PHONY: all

all: $(TARGET_bin) $(TARGET_STATIC_bin) $(TARGET_lib) $(TARGET_STATIC_lib)

$(TARGET_bin) : geo.obj test.obj simple.obj service.obj
	$(LINK_COMMAND_LINE) $^ -o $@

$(TARGET_STATIC_bin) : geo.obj test.obj simple_static.obj service.obj
	$(LINK_COMMAND_LINE) $^ -o $@

$(TARGET_lib): pyAPI.obj geo.obj simple.obj service.obj
	$(LINK_COMMAND_LINE) -shared $^ -o $@

$(TARGET_STATIC_lib): pyAPI_static.obj geo.obj simple_static.obj service.obj
	$(LINK_COMMAND_LINE) -shared $^ -o $@

simple.obj: src/simple.cu
	$(NVCC_COMMAND_LINE) $< -o $@

simple_static.obj: src/simple.cu
	$(NVCC_COMMAND_LINE) -DSTATIC $< -o $@

service.obj : src/service.cu
	$(NVCC_COMMAND_LINE) $< -o $@

geo.obj: src/geo.cpp
	$(GCC_COMMAND_LINE) $< -o $@

test.obj: src/test.cpp
	$(GCC_COMMAND_LINE) -fopenmp $< -o $@

pyAPI.obj: src/pyAPI.cpp
	$(GCC_PY_COMMAND_LINE) $< -o $@ -DPYMODNAME=pyEmpMe

pyAPI_static.obj: src/pyAPI.cpp
	$(GCC_PY_COMMAND_LINE) $< -o $@ -DSTATIC -DPYMODNAME=pyEmpMe_static

clean:
	rm -f geo.obj test.obj simple_static.obj service.obj simple.obj pyAPI.obj pyAPI_static.obj $(TARGET_STATIC_bin) $(TARGET_STATIC_lib) $(TARGET_bin) $(TARGET_lib)
