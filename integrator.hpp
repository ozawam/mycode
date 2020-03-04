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


void leap_frog(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          PS::F32 & time_sys,
          const PS::S64 nloop,
          PS::DomainInfo & dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) ;



void hermite(PS::ParticleSystem<FPGrav> & system_grav,
          PS::F64 & dt,
          PS::F32 & time_sys,
          const PS::S64 n_loop,
          PS::DomainInfo & dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
);

void initial_timestep(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys,PS::F32 & time_sys);





