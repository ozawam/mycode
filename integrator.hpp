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
/*
template<class Tpsys>
void kick(Tpsys & system ,const PS::F64 dt); 

template<class Tpsys>
void drift(Tpsys & system ,const PS::F64 dt); 
*/
/*
template<class Tpsys>
void kick(const Tpsys & ,const PS::F64 dt); 

template<class Tpsys>
void drift(const Tpsys & ,const PS::F64 dt); 
*/

template<class Tpsys>
void kick(Tpsys & system, 
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

