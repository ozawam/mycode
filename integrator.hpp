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

//template<class Tpsys>
void kick(PS::ParticleSystem<FPGrav> & system, 
           const PS::F64 dt) ;

//template<class Tpsys>
void drift(PS::ParticleSystem<FPGrav> & system, 
           const PS::F64 dt) ;

//#ifdef MULTI_WALK
void leap_frog(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          const PS::F32 time_sys,
          const PS::S64 nloop,
          PS::DomainInfo dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::S64 n_tot,
          const PS::F32 theta,
          const PS::S32 n_leaf_limit,
          const PS::S32 n_group_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) ;
/*
#else
void leap_frog(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          const PS::F32 time_sys,
          const PS::S64 nloop,
          PS::DomainInfo dinfo,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav
) ;*/
