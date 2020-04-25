#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include <time.h>                                                                                                                                        
#include <sys/time.h>
#include "tools.h"
#include "output.h"

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);

double E0;
struct reb_vec3d L0;
//double _rho_planet = 2.0;//g/cm^3
//double rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c);

//global const                                                                                                                                           
double M_sun = 1.989*pow(10.0,33.0);//g
double M_earth = 5.9724*pow(10.0,27.0);//g
double one_au=1.49597870*pow(10.0,13.0);//cm
double one_year = 365.0*24.0*60.0*60.0;//year-s

  
//global parameter
double tmax = 1.0e7 * 2.0 * M_PI;
double min_dt = 1.0e-22 * 2.0 * M_PI;
double C_aero_drag = 1.0;
double _rho_planet = 2.0;//g/cm^3

int main(int argc, char* argv[]){

  printf("Creating simulation from binary file and integrating until tmax\n");
  char PAST_SNAP[50];
  int restart_number = atoi(argv[1]) ;
  sprintf(PAST_SNAP,"../output/restart/restart%05d.bin",restart_number);
        struct reb_simulation* r = reb_create_simulation_from_binary(PAST_SNAP);
        // Need to reset function pointers
//  printf("check 1\n");

//Setup constants.
    r->dt                   = 0.01*2.0*M_PI;                // initial timestep


// Setup callback function for collision.
    r->collision_resolve = reb_collision_resolve_merge_pass_through;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;

// Setup integrator.
  //   r->integrator           = REB_INTEGRATOR_IAS15;
     r->integrator    = REB_INTEGRATOR_MERCURIUS;
     r->ri_mercurius.hillfac = 1;
  //   r->ri_ias15.epsilon         = 1e-10;
  //     r->ri_ias15.iterations_max_exceeded = 1000;
  //    r->ri_ias15.min_dt      = 1e-5 * 2.0 * M_PI;
     r->ri_ias15.min_dt      = min_dt;

// Setup callback function for velocity dependent forces.
   r->additional_forces     = additional_forces; 
   r->force_is_velocity_dependent = 1;

//  printf("check 2\n");
// Setup callback function for outputs.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = restart_number + 1;


        reb_integrate(r, tmax);        
        printf("Done.\n");
}


void additional_forces(struct reb_simulation* const r){
    // velocity dependent drag force.
	struct reb_particle* const particles = r->particles;
	//Setup constants. Units AU,Msun,T G=1
        double L_unit=one_au;
        double M_unit=M_sun;
        double G_unit=6.67408*pow(10.0,-8.0);
        double T_unit=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));

	double rho_planet = _rho_planet * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
	double mass_sun = 1.0;// particles[0].m;
	double Cd = C_aero_drag; // coeficient aero drag
	
	//Calculation of Aero Drag
	const int N = r->N; //Sun + planetesimals
	for (int i=2;i<N;i++){
        //variable
	double r_sq = particles[i].x*particles[i].x + particles[i].y*particles[i].y;
	double inv_r = 1.0 / sqrt(r_sq);
	double ev[3] = {-particles[i].y*inv_r, particles[i].x*inv_r, 0.0}; // unit vector of kepler velocity
	double vkep[3] = { sqrt(mass_sun * inv_r) * ev[0],sqrt(mass_sun * inv_r) * ev[1],sqrt(mass_sun * inv_r) * ev[2],}; // kepler velocity //0 x, 1 y, 2 z
	double vgas[3] = {0.995*vkep[0],0.995*vkep[1],0.995*vkep[2]};//{(1.0 - eta)*vkep[0],(1.0 - eta)*vkep[1],(1.0 - eta)*vkep[2]};
	double u[3] = {particles[i].vx - vgas[0],particles[i].vy - vgas[1],particles[i].vz - vgas[2]};
	double m_50 = 4.0/3.0 * pow(5000.0/one_au,3.0) *M_PI * rho_planet;//r=50m 2019/06/12
	double rate_50 = 8.958600e+25/M_sun/m_50;
	double rplanet = cbrt(3.0*particles[i].m/rate_50/(4.0*M_PI*rho_planet)) * one_au; //unit:cm

	double C1 = 11.1 /(rplanet * rplanet) /(24.0 * 60.0 * 60.0 /T_unit)  ;//* pow(one_au,2.0); //unit:day^-1 -> T_unit^-1

	double sys_acc_gd[3] = {(-C1*Cd*u[0] ),(-C1*Cd*u[1] ) ,(-C1*Cd*u[2])  };

	particles[i].ax += sys_acc_gd[0];
        particles[i].ay += sys_acc_gd[1];
        particles[i].az += sys_acc_gd[2];
    }
}


int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c){
    // This function passes the collision to the default merging routine.
    // If a merger occured, that routine will return a value other than 0.
    // This function then outputs some information about the merger.
//    int result = reb_collision_resolve_hardsphere(r,c);//reb_collision_resolve_merge(r, c);
        int result =reb_collision_resolve_merge(r, c);
    if (result!=0){
        printf("A merger occured! Particles involved: %d, %d.\n",c.p1,c.p2);
    }
    return result;
}

void heartbeat(struct reb_simulation* r){
      //const
        double M_sun=1.989*pow(10.0,33.0);//g
        double one_au=1.49597870*pow(10.0,13.0);//cm
        double one_year = 365.0*24.0*60.0*60.0;//year-s
        double L_unit=one_au;
        double M_unit=M_sun;
        double G_unit=6.67408*pow(10.0,-8.0);

        double t_16=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));
        double T_unit=t_16;


     //   int snap_n = (int)(r->t)/(100.0*2.0*M_PI);
        int snap_n = (r->hash_ctr);
        char SNAP[50];
        char SNAP_xyz[50];
        char SNAP_restart[50];
    if (reb_output_check(r, 100.0*2.0*M_PI)){  
        reb_output_timing(r, tmax);
           }
if((r->t)*(T_unit/one_year) < 10000.0){
    if (reb_output_check(r, 100.0*2.0*M_PI)){
       // Once per 100 years, output the relative energy error to a dat file
        FILE* f = fopen("../energy.dat","a");
        reb_integrator_synchronize(r);
        double E = reb_tools_energy(r);  
        double relE = fabs((E-E0)/E0);

        struct reb_vec3d L = reb_tools_angular_momentum(r); 
        struct reb_vec3d relL = {fabs((L.x-L0.x)/L0.x),fabs((L.y-L0.y)/L0.y) , fabs((L.z-L0.z)/L0.z)};    

        fprintf(f,"%e\t %e\t %e\t %e\t %e\n",(r->t)*(T_unit/one_year), relE, relL.x, relL.y, relL.z);
        fclose(f);

//  printf("check 2\n");
       sprintf(SNAP,"../output/snap01/snap%05d.dat",snap_n);
       sprintf(SNAP_xyz,"../output/xyz_snap01/snap%05d.dat",snap_n);
       sprintf(SNAP_restart,"../output/restart/restart%05d.bin",snap_n);

       reb_output_orbits(r,SNAP);
       reb_output_ascii(r,SNAP_xyz);
       //For restart
       reb_output_binary(r, SNAP_restart);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }
}
else{
    if (reb_output_check(r, 100.0*2.0*M_PI)){
       // Once per 100 years, output the relative energy error to a dat file
        FILE* f = fopen("../energy.dat","a");
        reb_integrator_synchronize(r);
        double E = reb_tools_energy(r);  
        double relE = fabs((E-E0)/E0);

        struct reb_vec3d L = reb_tools_angular_momentum(r); 
        struct reb_vec3d relL = {fabs((L.x-L0.x)/L0.x),fabs((L.y-L0.y)/L0.y) , fabs((L.z-L0.z)/L0.z)};    

        fprintf(f,"%e\t %e\t %e\t %e\t %e\n",(r->t)*(T_unit/one_year), relE, relL.x, relL.y, relL.z);
        fclose(f);

//  printf("check 2\n");
       sprintf(SNAP,"../output/snap01/snap%05d.dat",snap_n);
       sprintf(SNAP_xyz,"../output/xyz_snap01/snap%05d.dat",snap_n);
       sprintf(SNAP_restart,"../output/restart/restart%05d.bin",snap_n);

       reb_output_orbits(r,SNAP);
       reb_output_ascii(r,SNAP_xyz);
       //For restart
       reb_output_binary(r, SNAP_restart);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }

}
}

