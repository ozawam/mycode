#include<particle_simulator.hpp>
#include "integrator.hpp"


#define SECONDORDER 5e-1
#define THIRDORDER  1.666666666666667e-1
#define FOURTHORDER 4.166666666666667e-2
#define FIFTHORDER  8.333333333333333e-3


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
          PS::F32 & time_sys,
          const PS::S64 n_loop,
          PS::DomainInfo & dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) {

    
        kick(system_grav, dt * 0.5);
        
        time_sys += dt;
        drift(system_grav, dt);
        
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        
//    std::cerr << "check3 " << std::endl;
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

        tree_grav.calcForceAllAndWriteBack(CalcGravity<FPGrav>,
                                           CalcGravity<PS::SPJMonopole>,
                                           system_grav,
                                           dinfo);

#endif
        
        kick(system_grav, dt * 0.5);
}

/*
//#######################################################
// Hermite
//#######################################################  

//1.predictor
// ---------------------------------------------------------------------------

static void predict(REAL t_gl, int ni, int *index, struct Particle *particle, struct Address *address)
{
  int i, iadr, padr, id;
  REAL dt, dt2, dt3;

  for(i = 0; i < ni; i++){
    iadr = index[i];
    id   = particle[iadr].id;
    dt   = t_gl - particle[iadr].t;
    dt2 = dt * dt;
    dt3 = dt * dt2;

    posvel[i].xpos = particle[iadr].xpos + particle[iadr].xvel * dt + particle[iadr].xacc * dt2 * SECONDORDER + particle[iadr].xjrk * dt3 * THIRDORDER;
    posvel[i].ypos = particle[iadr].ypos + particle[iadr].yvel * dt + particle[iadr].yacc * dt2 * SECONDORDER + particle[iadr].yjrk * dt3 * THIRDORDER;
    posvel[i].zpos = particle[iadr].zpos + particle[iadr].zvel * dt + particle[iadr].zacc * dt2 * SECONDORDER + particle[iadr].zjrk * dt3 * THIRDORDER;
    posvel[i].xvel = particle[iadr].xvel + particle[iadr].xacc * dt + particle[iadr].xjrk * dt2 * SECONDORDER;
    posvel[i].yvel = particle[iadr].yvel + particle[iadr].yacc * dt + particle[iadr].yjrk * dt2 * SECONDORDER;
    posvel[i].zvel = particle[iadr].zvel + particle[iadr].zacc * dt + particle[iadr].zjrk * dt2 * SECONDORDER;
    posvel[i].id   = (float)id;
  }
  return;
}





//2.a1j & jerk1j
// --------------------------------------------------------------------------- 





//3.corrector
// --------------------------------------------------------------------------- 




//4.next step
// --------------------------------------------------------------------------- 






// --------------------------------------------------------------------------- 
//  time integration
// --------------------------------------------------------------------------- 


*/
