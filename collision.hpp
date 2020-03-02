#pragma once
#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif

#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_gpu_cuda.hpp"
#endif

#include "user-defined.hpp"
#include "integrator.hpp"


void collision(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys,PS::F32 & time_sys, const PS::F64 rij);

