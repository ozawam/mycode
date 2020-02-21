#include<particle_simulator.hpp>
#include "integrator.hpp"


//#######################################################
// leap-frog
//#######################################################  

//template<class Tpsys>
void kick(PS::ParticleSystem<FPGrav> & system, 
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].vel  += system[i].acc * dt;
    }
}

//template<class Tpsys>
void drift(PS::ParticleSystem<FPGrav> & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos  += system[i].vel * dt;
    }
}

void leap_frog(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          PS::F32 time_sys,
          const PS::S64 n_loop,
          PS::DomainInfo dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::S64 n_tot,
          const PS::F32 theta,
          const PS::S32 n_leaf_limit,
          const PS::S32 n_group_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) {

    
//    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
//    if(time_sys == 0.0){
    std::cerr << "check1 " << std::endl;
//    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
//    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
//    std::cerr << "check2 " << std::endl;
//}
        kick(system_grav, dt * 0.5);
        
        time_sys += dt;
        drift(system_grav, dt);
        
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        
        system_grav.exchangeParticle(dinfo);

#ifdef MULTI_WALK
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                    RetrieveKernel,
                                                    tag_max,
                                                    system_grav,
                                                    dinfo,
                                                    n_walk_limit,
                                                    true);
#else
    std::cerr << "check2 " << std::endl;
        tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           system_grav,
                                           dinfo);

#endif
        
    std::cerr << "check3 " << std::endl;
        kick(system_grav, dt * 0.5);
}


//#######################################################
// Hermite
//#######################################################  













