#include <limits>
#include<particle_simulator.hpp>
#include "integrator.hpp"


#define SECONDORDER 5e-1
#define THIRDORDER  1.666666666666667e-1
#define FOURTHORDER 4.166666666666667e-2
#define FIFTHORDER  8.333333333333333e-3

//global const                                                                                                                                           
PS::F64 M_sun = 1.989*pow(10.0,33.0);//g
PS::F64 M_earth = 5.9724*pow(10.0,27.0);//g
PS::F64 one_au=1.49597870*pow(10.0,13.0);//cm
PS::F64 one_year = 365.0*24.0*60.0*60.0;//year-s

PS::F64 eta = 5.0e-2 ;
PS::F64 G = 1.0;

PS::S64 col_n;

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

static void predict(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys,PS::F32 & time_sys)
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

static void correct(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys, PS::F32 & time_sys)
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

static void timestep(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys,PS::F32 & time_sys)
{
  PS::S32 i;
  PS::S32 n = particle.getNumberOfParticleLocal();
//  PS::F64 eta = 3.0e-3 ;

  PS::F64vec a1;
  PS::F64vec j1;
  PS::F64vec s1;
  PS::F64vec c1;

  PS::F64 A1;
  PS::F64 J1;
  PS::F64 S1;
  PS::F64 C1;

  PS::F64 min_dt;

  for(i = 0; i < n; i++){

    a1 = particle[i].acc ;
    j1 = particle[i].jrk ;

    s1 = particle[i].s0 + ( particle[i].c0 * dt_sys ) ;
    c1 = particle[i].c0 ;
    
    A1 = sqrt(a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2] ) ;
    J1 = sqrt(j1[0] * j1[0] + j1[1] * j1[1] + j1[2] * j1[2] ) ;
    S1 = sqrt(s1[0] * s1[0] + s1[1] * s1[1] + s1[2] * s1[2] ) ;
    C1 = sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2] ) ;

    particle[i].pdt = eta * sqrt((A1 * S1 + J1 * J1 )/(J1 * C1 + S1 * S1 ))  ;
  
    if(i == 0){ 
      min_dt = particle[i].pdt  ;   
      }
    else if( particle[i].pdt < min_dt  ){
      min_dt = particle[i].pdt  ;  
      } 
//finding Sun for next output
// --------------------------------------------------------------------------- 
    if(particle[i].mass > 0.9){
      particle[0].id_sun = i;
      }
  }

//  dt_sys = min_dt;
      dt_sys = PS::Comm::getMinValue(min_dt);

  return;
}


void initial_timestep(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys,PS::F32 & time_sys)
{
  PS::S32 i;
  PS::S32 n = particle.getNumberOfParticleLocal();
//  PS::F64 eta = 3.0e-3 ;

  PS::F64vec a1;
  PS::F64vec j1;

  PS::F64 A1;
  PS::F64 J1;

  PS::F64 min_dt;

  for(i = 0; i < n; i++){

    a1 = particle[i].acc ;
    j1 = particle[i].jrk ;
 
    A1 = sqrt(a1[0] * a1[0] + a1[1] * a1[1] + a1[2] * a1[2] ) ;
    J1 = sqrt(j1[0] * j1[0] + j1[1] * j1[1] + j1[2] * j1[2] ) ;

    particle[i].pdt = eta * A1/J1  ;
//   fprintf(stdout, "timestep: %10.7f  A1:%10.7f  J1:%10.7f  j1[0]:%10.7f   \n", particle[i].pdt , A1, J1, j1[0] );
    if(i == 0){ 
      min_dt = particle[i].pdt  ;   
      }
    else if( particle[i].pdt < min_dt  ){
      min_dt = particle[i].pdt  ;  
      } 
  //    fprintf(stdout, "timestep: %10.7f \n", min_dt);
//finding Sun for next output
// --------------------------------------------------------------------------- 
    if(particle[i].mass > 0.9){
      particle[0].id_sun = i;
      }

  }

    if(min_dt != std::numeric_limits<double>::infinity() ){ 
//      dt_sys = min_dt;
      dt_sys = PS::Comm::getMinValue(min_dt);
        }
  return;
}



//#######################################################
// Merge
//#######################################################  

//Merge
// ---------------------------------------------------------------------------

void merge(PS::ParticleSystem<FPGrav> & particle,  const PS::F32 time_sys)
{

    PS::S32 n_remove = 0;

    col_n = particle[0].collisions_N;
    std::vector<PS::S32> idx(col_n);
//    PS::S32 idx[col_n];
  for(PS::S64 ci = 0; ci < particle[0].collisions_N ; ci++){
//     fprintf(stdout, "CHECK1 \n");

    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    PS::S64 i = particle[0].COL_P[2*ci];
    PS::S64 j = particle[0].COL_P[2*ci+1];   //want j to be removed particle
//    if (j<i){
//        i = particle[0].COL_P[2*ci+1];
//        j = particle[0].COL_P[2*ci];
//    }
    fprintf(stdout, "j: %lld \n", j);
//    idx.push_back(j);
    idx[ci] = j ;
    n_remove ++;

    // Check collision of j particle
    for(PS::S64 iC = ci+1 ; iC < particle[0].collisions_N ; iC++){
    if(particle[0].COL_P[2*iC] == j  ){
      particle[0].COL_P[2*iC] = i ;
    }
    else if(particle[0].COL_P[2*iC+1] == j  ){
      particle[0].COL_P[2*iC+1] = i ;
    }

    }


    // Scale out energy from collision - initial energy
    PS::F64 energy_offset = 0.0;
    PS::F64 Ei=0.0, Ef=0.0;
    
        {
            // Calculate energy difference in inertial frame
            Ei = 0.5 * particle[i].mass * ( particle[i].vel * particle[i].vel ) + 0.5 * particle[j].mass * ( particle[j].vel * particle[j].vel ) ;
        }
        
    //merge or rebound
   	//constant
        PS::F64 L_unit=one_au;
        PS::F64 M_unit=M_sun;
        PS::F64 G_unit=6.67408*pow(10.0,-8.0);

        PS::F64 t_16=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));
        PS::F64 T_unit=t_16;

        PS::F64vec vr;
        vr = particle[j].vel - particle[i].vel ;

        PS::F64 eps = 0.7;
        PS::F64 vreb;
        vreb = fabs( eps * sqrt( vr * vr ));

        PS::F64  vesc12;
        PS::F64  r1,r2;
        PS::F64  _rho_planet = 2.0;//g/cm^3
        PS::F64  rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
        r1 = pow(3.0*particle[i].mass/(4.0*M_PI*rho_planet),1.0/3.0);
        r2 = pow(3.0*particle[j].mass/(4.0*M_PI*rho_planet),1.0/3.0);
        vesc12 = sqrt(2.0*G*(particle[i].mass + particle[j].mass)/(r1 + r2));

        char filename[1024] = "result/collision_after.dat";
        FILE* of = fopen(filename,"a"); //change mode w->a 
        if (of==NULL){
        fprintf(stdout,"Can not open file.");
        }

        fprintf(of,"%le\t%e\t%e\t%e\t\n",time_sys,vreb,vesc12,Ei);
        //time T_unit

        particle[i].writeAscii(of);
        fprintf(of,"%lld\t%e\t\n",i,(particle[i].mass)*M_unit);
        particle[j].writeAscii(of);
        fprintf(of,"%lld\t%e\t\n",j,(particle[j].mass)*M_unit);

        // Merge by conserving mass, volume and momentum
        particle[i].pos = ((particle[i].pos * particle[i].mass) + (particle[j].pos * particle[j].mass))/(particle[i].mass + particle[j].mass) ;
        particle[i].vel = ((particle[i].vel * particle[i].mass) + (particle[j].vel * particle[j].mass))/(particle[i].mass + particle[j].mass) ;
        particle[i].mass  = particle[i].mass + particle[j].mass ;
        particle[i].id = i;

        particle[i].writeAscii(of);
        fprintf(of,"%lld\t%e\t%e\t\n",i,(particle[i].mass)*M_unit,(time_sys)*(T_unit/one_year));        
        //writeAscii: time T_unit~0.15925 year, mass Msun, length AU
        //footer: time year, mass g :::::conversion
   
        if(vreb < vesc12){
        fprintf(of,"merge\n");
        fclose(of);  
        }
        else{
        fprintf(of,"rebound\n"); 
        fclose(of); 
        }

        // Keeping track of energy offst
        {  
            Ef = 0.5 * particle[i].mass * ( particle[i].vel * particle[i].vel ) ;
        }
        energy_offset += Ei - Ef;

    }

        // Removing particle j
//        std::sort(idx.begin(),idx.end(),std::greater<PS::S32>());
//        particle.removeParticle( idx, n_remove);
           //Keeping sort when removing particle 2020/03/02   
		        PS::S32 N = particle.getNumberOfParticleGlobal();
            std::sort(idx.begin(),idx.end(),std::greater<PS::S32>());
            for(PS::S64 I=0; I<n_remove; I++){
		        for(PS::S64 J=idx[I]; J<N-1; J++){
			      particle[J] = particle[J+1];
		        particle[J].id = J ;
            }
            N--;
            particle[0].N_active --;

            particle[N].mass = 0.0 ;
            particle[N].pos[2] =1.0e+23;

//            PS::S32 id_last[1] = {N-1};
//            particle.removeParticle( id_last, 1);
            }

//            PS::S32 id_last[1] = {71};
//            particle.removeParticle( id_last, 1);
            particle[0].collisions_N =0;
            particle[0].COL_P.erase( particle[0].COL_P.begin() ,  particle[0].COL_P.end() );
  return;
}


//#######################################################
// Gas drag
//#######################################################  

static void gas_drag(PS::ParticleSystem<FPGrav> & particle,  PS::F64 & dt_sys, PS::F32 & time_sys)
{
//gas
// ---------------------------------------------------------------------------
	//Setup constants related units. Units AU,Msun,T G=1
  PS::F64 L_unit=one_au;
  PS::F64 M_unit=M_sun;
  PS::F64 G_unit=6.67408*pow(10.0,-8.0);
  PS::F64 T_unit=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));
//const                                                                                                                                           
  PS::F64 C_aero_drag = 1.0;
  PS::F64 _rho_planet = 2.0;//g/cm^3
	PS::F64 rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
	PS::F64 mass_sun = 1.0;// particles[0].m;
	PS::F64 Cd = C_aero_drag; // coeficient aero drag
	PS::F64 alpha = 0.995;

  PS::S32 i;
  PS::S32 n = particle.getNumberOfParticleLocal();
  
  for(i = 0; i < n; i++){
if(particle[i].mass*M_sun <1.0E+30 ){
	//Calculation of Aero Drag
/*//variable
  PS::F64 gx = particle[i].pos[0];
  PS::F64 gy = particle[i].pos[1];
  PS::F64 gz = particle[i].pos[2];
  PS::F64 gm = particle[i].mass;
  PS::F64vec xi = particle[i].pos;
  PS::F64vec vi = particle[i].vel;
  PS::F64vec ai = particle[i].acc;
*/

//variable
  PS::F64vec pos_s =particle[particle[0].id_sun].pos;
  PS::F64vec vel_s =particle[particle[0].id_sun].vel;
  PS::F64vec acc_s =particle[particle[0].id_sun].vel;

  PS::F64vec xi = particle[i].pos - pos_s;
  PS::F64vec vi = particle[i].vel - vel_s;
  PS::F64vec ai = particle[i].acc - acc_s;
  PS::F64 gx = xi[0];
  PS::F64 gy = xi[1];
  PS::F64 gz = xi[2];
  PS::F64 gm = particle[i].mass;

	PS::F64 r_sq = gx*gx + gy*gy;
	PS::F64 inv_r = 1.0 / sqrt(r_sq);
	PS::F64vec ev = {-gy*inv_r, gx*inv_r, 0.0}; // unit vector of kepler velocity 
	PS::F64vec vkep = sqrt(mass_sun * inv_r) * ev; // kepler velocity //0 x, 1 y, 2 z
	PS::F64vec vgas = alpha*vkep;//{(1.0 - eta)*vkep[0],(1.0 - eta)*vkep[1],(1.0 - eta)*vkep[2]};
	PS::F64vec u = vi - vgas;


  PS::F64 Cg = alpha*sqrt(G_unit*mass_sun);
  PS::F64 r_deriv = inv_r * (gx * vi[0] + gy * vi[1])  ;
	PS::F64vec vkep_deriv = {(-3.0/2.0)*pow(inv_r,5.0/2.0)*r_deriv*(-1)*gy + pow(inv_r,3.0/2.0)*(-1)*vi[1] , (-3.0/2.0)*pow(inv_r,5.0/2.0)*r_deriv*gx + pow(inv_r,3.0/2.0)*vi[0]   , 0 };
  PS::F64vec u_deriv = ai - Cg*vkep_deriv;


//radius and C1  
//	PS::F64 m_50 = 4.0/3.0 * pow(5000.0/one_au,3.0) *M_PI * rho_planet;//r=50m 2019/06/12
//	PS::F64 rate_50 = 8.958600e+25/M_sun/m_50;
//	PS::F64 rplanet = cbrt(3.0*gm/rate_50/(4.0*M_PI*rho_planet)) * one_au; //unit:cm(Unit is cm in calculation of C1)
	PS::F64 rplanet = 10.0 * 100.0; //unit:cm(Unit is cm in calculation of C1)

	PS::F64 C1 = 11.1 /(rplanet * rplanet) /(24.0 * 60.0 * 60.0 /T_unit)  ;//* pow(one_au,2.0); //unit:day^-1 -> T_unit^-1, unit:cm -> AU  


	PS::F64vec sys_acc_gd = (-C1*Cd*u );
	PS::F64vec sys_jrk_gd = (-C1*Cd*u_deriv );

        particle[i].acc += sys_acc_gd;
        particle[i].jrk += sys_jrk_gd;
}//end of if
}//end of for
}


// --------------------------------------------------------------------------- 
//  time integration
// --------------------------------------------------------------------------- 

void hermite(PS::ParticleSystem<FPGrav> & system_grav,
          PS::F64 & dt,
          PS::F32 & time_sys,
          const PS::S64 n_loop,
          PS::DomainInfo & dinfo,
          const PS::S32 tag_max,
          const PS::S32 n_walk_limit,
          PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole & tree_grav
) {

//1.predictor
// ---------------------------------------------------------------------------
            system_grav[0].hermite_step =1;
  predict(system_grav, dt , time_sys);

//2.a1j & jerk1j
// --------------------------------------------------------------------------- 
//            system_grav[0].hermite_step =2;
            system_grav[0].hermite_step =0;//no collision
        
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        
//    std::cerr << "check3 " << std::endl;
        system_grav.exchangeParticle(dinfo);

  tree_grav.calcForceAllAndWriteBack(CalcGravity(),
                                     system_grav,
                                     dinfo);

  gas_drag(system_grav, dt , time_sys);

//!collision 
// ---------------------------------------------------------------------------

//3.1 corrector
// --------------------------------------------------------------------------- 
  correct(system_grav, dt , time_sys);


//!merge
// ---------------------------------------------------------------------------
//  merge(system_grav, time_sys); 
if(system_grav[0].collisions_N > 0){  
  if(PS::Comm::getRank() == 0) {
  merge(system_grav, time_sys); 

  fprintf(stdout, "CHECK_COL_end \n");
  }  
   else {
  system_grav.setNumberOfParticleLocal(0);
  fprintf(stdout, "CHECK_else \n");
    }
}


//3.2 new a1j & jerk1j 
// --------------------------------------------------------------------------- 
            system_grav[0].hermite_step =3;
        system_grav.exchangeParticle(dinfo);
        tree_grav.calcForceAllAndWriteBack(CalcGravity(),
                                           system_grav,
                                           dinfo);
  gas_drag(system_grav, dt , time_sys);

//4.next step
// --------------------------------------------------------------------------- 
            system_grav[0].hermite_step =4;
        time_sys += dt;
if(col_n == 0){
  timestep(system_grav, dt, time_sys);
}

else{
    initial_timestep( system_grav, dt, time_sys);
}

}

