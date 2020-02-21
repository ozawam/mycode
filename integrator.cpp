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

void leap_frog1(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          PS::F32 time_sys,
          const PS::S64 n_loop
//          PS::DomainInfo dinfo
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
        
//        if(n_loop % 4 == 0){
//            dinfo.decomposeDomainAll(system_grav);
//        }
        
//        system_grav.exchangeParticle(dinfo);
}
        

void leap_frog2(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          PS::F32 time_sys,
          const PS::S64 n_loop
//          PS::DomainInfo dinfo
) {
        kick(system_grav, dt * 0.5);
}


//#######################################################
// Hermite
//#######################################################  













