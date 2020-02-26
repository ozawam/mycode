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

        tree_grav.calcForceAllAndWriteBack(CalcGravity(),
//                                           CalcGravity(),
                                           system_grav,
                                           dinfo);

#endif
        
        kick(system_grav, dt * 0.5);
}


//#######################################################
// Hermite
//#######################################################  

//1.predictor
// ---------------------------------------------------------------------------

static void predict(PS::ParticleSystem<FPGrav> & particle, const PS::F64 & dt_sys,PS::F32 & time_sys)
{
  PS::S32 i,iadr,id;
  PS::F64 dt, dt2, dt3;
  PS::S32 n = particle.getNumberOfParticleLocal();

  for(i = 0; i < n; i++){
    iadr = particle[i].id;
    id   = particle[iadr].id;
    dt   = dt_sys;
    dt2 = dt * dt;
    dt3 = dt * dt2;

    particle[i].x0 = particle[i].pos ;
    particle[i].v0 = particle[i].vel ;
    particle[i].a0 = particle[i].acc ;
    particle[i].j0 = particle[i].jrk ;


    particle[i].xp = particle[i].x0 + particle[i].v0 * dt + particle[i].a0 * dt2 * SECONDORDER + particle[i].j0 * dt3 * THIRDORDER;

    particle[i].vp = particle[i].v0 + particle[i].a0 * dt + particle[i].j0 * dt2 * SECONDORDER;
    
    particle[i].pos =  particle[i].xp ;
    particle[i].vel =  particle[i].vp ;
//    posvel[i].id   = (double)id;
  }
  return;
}




//2.a1j & jerk1j -> in void hermite
// --------------------------------------------------------------------------- 


//3.1 corrector
// --------------------------------------------------------------------------- 

static void correct(PS::ParticleSystem<FPGrav> & particle, const PS::F64 & dt_sys, PS::F32 & time_sys)
{

  PS::F64vec xc;
  PS::F64vec vc;

  PS::S32 i,iadr,id;
  PS::S32 n = particle.getNumberOfParticleLocal();
  PS::F64 dt, dtsinv;

  dt    = dt_sys;
  dtsinv = 1.0 / dt;

  for(i = 0; i < n; i++){
    particle[i].a1 = particle[i].acc ;
    particle[i].j1 = particle[i].jrk ;

    particle[i].s0  = 2.0 * (- 3.0 * (particle[i].a0 - particle[i].a1) - (2.0 * particle[i].j0 + particle[i].j1) * dt) * dtsinv * dtsinv;
    particle[i].c0  = 6.0 * (2.0 * (particle[i].a0 - particle[i].a1) + (particle[i].j0 + particle[i].j1) * dt) * dtsinv * dtsinv * dtsinv;
    xc = particle[i].xp + particle[i].s0 * dt * dt * dt * dt * FOURTHORDER + particle[i].c0 * dt * dt * dt * dt * dt * FIFTHORDER;
    vc = particle[i].vp + particle[i].s0 * dt * dt * dt * THIRDORDER + particle[i].c0 * dt * dt * dt * dt * FOURTHORDER;
  
    particle[i].pos =  xc ;
    particle[i].vel =  vc ;
  }

  return;
}

//3.2 new a1j & jerk1j -> in void hermite
// --------------------------------------------------------------------------- 

//4.next step
// --------------------------------------------------------------------------- 

static void timestep(PS::ParticleSystem<FPGrav> & particle, const PS::F64 & dt_sys,PS::F32 & time_sys)
{
  PS::S32 i,iadr,id;
  PS::F64 dt, dt2, dt3;
  PS::S32 n = particle.getNumberOfParticleLocal();

  for(i = 0; i < n; i++){
    iadr = particle[i].id;
    id   = particle[iadr].id;
    dt   = dt_sys;
    dt2 = dt * dt;
    dt3 = dt * dt2;

    particle[i].x0 = particle[i].pos ;
    particle[i].v0 = particle[i].vel ;
    particle[i].a0 = particle[i].acc ;
    particle[i].j0 = particle[i].jrk ;


    particle[i].xp = particle[i].x0 + particle[i].v0 * dt + particle[i].a0 * dt2 * SECONDORDER + particle[i].j0 * dt3 * THIRDORDER;

    particle[i].vp = particle[i].v0 + particle[i].a0 * dt + particle[i].j0 * dt2 * SECONDORDER;
    
    particle[i].pos =  particle[i].xp ;
    particle[i].vel =  particle[i].vp ;
//    posvel[i].id   = (double)id;
  }
  return;
}






// --------------------------------------------------------------------------- 
//  time integration
// --------------------------------------------------------------------------- 

void hermite(PS::ParticleSystem<FPGrav> & system_grav,
          const PS::F64 dt,
          PS::F32 & time_sys,
          const PS::S64 n_loop,
          PS::DomainInfo & dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) {

//1.predictor
// ---------------------------------------------------------------------------
  predict(system_grav, dt , time_sys);

//2.a1j & jerk1j
// --------------------------------------------------------------------------- 
        
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        
//    std::cerr << "check3 " << std::endl;
        system_grav.exchangeParticle(dinfo);

  tree_grav.calcForceAllAndWriteBack(CalcGravity(),
                                     system_grav,
                                     dinfo);

//3.1 corrector
// --------------------------------------------------------------------------- 
  correct(system_grav, dt , time_sys);

//3.2 new a1j & jerk1j 
// --------------------------------------------------------------------------- 
        system_grav.exchangeParticle(dinfo);
  tree_grav.calcForceAllAndWriteBack(CalcGravity(),
                                     system_grav,
                                     dinfo);

//4.next step
// --------------------------------------------------------------------------- 


        time_sys += dt;

}

