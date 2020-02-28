PS_PATH = ../../../FDPS/src/
INC = -I$(PS_PATH)
EXEFILE = nbody.out  

#CC = time g++
CC = g++
#CC = time mpig++
CFLAGS = -std=c++11 -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#use_phantom_grape_x86 = yes
#use_gpu_cuda = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

all:$(EXEFILE)

integrator.o:integrator.cpp integrator.hpp
	$(CC) $(INC) $(CFLAGS) -c $<  -o $@
OBJS = integrator.o integrator.hpp

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(PS_PATH)/phantom_grape_x86/G5/newton/libpg5
INC += -I$(PG_ROOT)
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   = -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_BUILD =
PG_CLEAN = 
endif

ifeq ($(use_gpu_cuda),yes)
CUDA_HOME = /usr/local/cuda
#CUDA_HOME = /gwfefs/opt/x86_64/cuda/7.5
NVCC = time $(CUDA_HOME)/bin/nvcc -std=c++11 -Xcompiler="-std=c++11 -O3"
INC  += -I$(CUDA_HOME)/samples/common/inc/
CFLAGS  += -DENABLE_GPU_CUDA
CLIBS = -L$(CUDA_HOME)/lib64 -lcudart -lgomp
force_gpu_cuda.o:force_gpu_cuda.cu
	$(NVCC) $(INC) -c -o $@ $<
OBJS = force_gpu_cuda.o
endif

nbody.out:nbody.cpp $(OBJS)
	$(PG_BUILD)
	$(CC) $(INC) $(CFLAGS) -o $@ $^ $(CLIBS)


clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf $(EXEFILE) *.o


distclean: clean
	$(PG_CLEAN)
	rm -f nbody.out
	rm -rf result

test: 
	# This command is only for FDPS developers.
	./test.py

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
