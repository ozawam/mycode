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

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c);

//global const                                                                                                                                           
double M_sun = 1.989*pow(10.0,33.0);//g
double M_earth = 5.9724*pow(10.0,27.0);//g
double one_au=1.49597870*pow(10.0,13.0);//cm
double one_year = 365.0*24.0*60.0*60.0;//year-s

  
//global parameter
double tmax = 1.0e-16 * 2.0 * M_PI;
double min_dt = 1.0e-15 * 2.0 * M_PI;
double C_aero_drag = 1.0;
double _rho_planet = 2.0;//g/cm^3



int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();

//Setup constants.
    r->dt                   = 0.01*2.0*M_PI;                // initial timestep
    r->G         = 1;

//Setup boundary.
   r->boundary    = REB_BOUNDARY_OPEN;
//   r->boundary    = REB_BOUNDARY_NONE;

// Setup integrator.
   r->integrator           = REB_INTEGRATOR_IAS15;
//   r->integrator    = REB_INTEGRATOR_MERCURIUS;
//   r->ri_mercurius.hillfac = 1;
//   r->ri_ias15.epsilon         = 1e-10;
//  r->ri_ias15.iterations_max_exceeded = 1000;
//  r->ri_ias15.min_dt      = min_dt;
//   r->integrator    = REB_INTEGRATOR_SABA;
//   r->integrator           = REB_INTEGRATOR_LEAPFROG;
//   r->integrator    = REB_INTEGRATOR_WHFAST;
//   r->ri_whfast.corrector = 11;    // Turn on symplectic correctors (11th order).
//   r->ri_whfast.safe_mode = 0;     // Turn of safe mode. Need to call integrator_synchronize() before outputs.

//Setup gravity solvers.
//    r->gravity    = REB_GRAVITY_TREE;
//    r->gravity    = REB_GRAVITY_COMPENSATED;
//    r->opening_angle2    = 0.5;        // This constant determines the accuracy of the tree code gravity estimate.
    const double boxsize = 1000.0;
    reb_configure_box(r,boxsize,2,2,1);

    // Initialize MPI
     // This can only be done after reb_configure_box.
//      reb_mpi_init(r);

// Setup collision.
//    r->collision            = REB_COLLISION_NONE;
    r->collision            = REB_COLLISION_DIRECT;
//    r->collision            = REB_COLLISION_LINE;
//    r->collision    = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_merge_pass_through;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
//    r->collision    = REB_COLLISION_TREE;
//    r->collision_resolve    = reb_collision_resolve_hardsphere;        // Choose merger collision routine.
//    r->collision_resolve    = reb_collision_resolve_merge;  
//   r->softening     = 0.02;        // Gravitational softening length

// Setup callback function for velocity dependent forces.
//   r->additional_forces     = additional_forces; 
//   r->force_is_velocity_dependent = 1;

// Setup callback function for outputs.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = 0;
//    r->usleep        = 10000;        // Slow down integration (for visualization only)

 //   double boxsize = 3.0;
 //   reb_configure_box(r,boxsize,1,1,1);


   // Star
    struct reb_particle star = {0};
    star.m = 1.0;
    star.r = 0.00465;
    //star.r = 0.1;
//if (r->mpi_id==0){
    reb_add(r, star);
//}

   // Planetesimal disk parameters
	int N_planetesimals = 1;    
	double planetesimal_mass = 1.0e-30;//total_disk_mass/N_planetesimals;//M_sun //unit:Msun;
    // Generate Planetesimal Disk
     for (int i=1;i<= N_planetesimals;i++){
        struct reb_particle pt = {0};
        double a  = 1.0;
        double e    = 0.0;// reb_random_rayleigh(0.003972877/sqrt(2.0));//circular orbit
        double inc  = 0.001;//reb_random_rayleigh(0.003972877/2.0/sqrt(2.0));//two-dimention
        double Omega = reb_random_uniform(0,2.*M_PI);
//        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi     = reb_random_uniform(0,2.*M_PI);

        double rho = _rho_planet;
  
        double radius_f = 1.0;//fold enlargement(radius enhancement)

        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis, phi); //change 3d to 2d in 2019/05/20
        pt.m = planetesimal_mass;
        pt.r = pow((3.0*planetesimal_mass*M_sun)/(4.0*M_PI*rho) ,1.0/3.0)/ one_au * radius_f; //unit is AU.
        pt.lastcollision = 0;
        reb_add(r, pt);
    }
//    r->N_active = r->N;
//    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
//    E0 = reb_tools_energy(r);
//    L0 = reb_tools_angular_momentum(r);

//    system("rm -rf energy.dat");
//    system("rm -rf output");
//    system("rm -rf collision_after.dat collision_before.dat");
//    system("mkdir output");
//    system("mkdir output/snap01");
//    system("mkdir output/xyz_snap01");

  //  reb_integrate(r, INFINITY);
    reb_integrate(r, tmax);

   // Cleanup
//    reb_mpi_finalize(r);
    reb_free_simulation(r);
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
	double rplanet = cbrt(3.0*particles[i].m/rate_50/(4.0*M_PI*rho_planet)) * one_au; //unit:cm(Unit is cm in calculation of C1)

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
//    int result = reb_collision_resolve_hardsphere(r,c);
	int result =reb_collision_resolve_merge(r, c);
    if (result!=0){
        printf("A merger occured! Particles involved: %d, %d.\n",c.p1,c.p2);
    }
    return result;
}

void heartbeat(struct reb_simulation* r){
      //const
//        double M_sun=1.989*pow(10.0,33.0);//g
//        double one_au=1.49597870*pow(10.0,13.0);//cm
//        double one_year = 365.0*24.0*60.0*60.0;//year-s
        double L_unit=one_au;
        double M_unit=M_sun;
        double G_unit=6.67408*pow(10.0,-8.0);

        double T_unit=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));

     //   int snap_n = (int)(r->t)/(100.0*2.0*M_PI);
        int snap_n = (r->hash_ctr);
//        char SNAP[50];
        char SNAP_xyz[50];
//        char SNAP_restart[50];

    if (reb_output_check(r, 10000.0*2.0*M_PI)){  
        reb_output_timing(r, tmax);
        //relative energy error 
//        double E = reb_tools_energy(r);
//        double relE = fabs((E-E0)/E0);

        //get orbital elements
//        struct reb_particle p = r->particles[1];
//        struct reb_particle star = r->particles[0];
//        struct reb_orbit o = reb_tools_particle_to_orbit(r->G,p,star);
//        printf("a2=%f,dE=%e,N=%d\n",o.a,relE,r->N);
           }

    if (reb_output_check(r, 1000000000.0*2.0*M_PI)){
       // Once per 100 years, output the relative energy error to a dat file
//        FILE* f = fopen("energy.dat","a");
//        reb_integrator_synchronize(r);
//        double E = reb_tools_energy(r);  
//        double relE = fabs((E-E0)/E0);

//        struct reb_vec3d L = reb_tools_angular_momentum(r); 
//        struct reb_vec3d relL = {fabs((L.x-L0.x)/L0.x),fabs((L.y-L0.y)/L0.y) , fabs((L.z-L0.z)/L0.z)};    

//        fprintf(f,"%e\t %e\t %e\t %e\t %e\n",(r->t)*(T_unit/one_year), relE, relL.x, relL.y, relL.z);
//        fclose(f);

//       sprintf(SNAP,"output/snap01/snap%05d.dat",snap_n);
       sprintf(SNAP_xyz,"output/snap%05d.dat",snap_n);
//       sprintf(SNAP_restart,"output/restart/restart%05d.bin",snap_n);  

       //For restart
//       reb_output_binary(r, SNAP_restart);

//       reb_output_orbits(r,SNAP);
       reb_output_ascii(r,SNAP_xyz);
   //    reb_output_orbits(r,"output/snap01/snap.dat");
       r->hash_ctr ++;
    }
}

