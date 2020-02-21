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
double ss_pos[3][3] =
    {
     {1.0, 3.0, 0.0},
     {-2.0, -1.0, 0.0}, 
     {1.0, -1.0, 0.0} 
};
double ss_vel[3][3] =
    {
     {0.0, 0.0, 0.0},
     {0.0, 0.0, 0.0},
     {0.0, 0.0, 0.0}
};

double ss_mass[6] =
    {
     3.0,
     4.0,
     5.0 
};                                                                                                                                           
  
//global parameter
double tmax = 100.0;



int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();

//Setup constants.
    r->dt                   = 0.00001;                // initial timestep
    r->G         = 1;

//Setup boundary.
//   r->boundary    = REB_BOUNDARY_OPEN;
   r->boundary    = REB_BOUNDARY_NONE;

// Setup integrator.
   r->integrator           = REB_INTEGRATOR_IAS15;

// Setup collision.
//    r->collision            = REB_COLLISION_NONE;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge_pass_through;
//    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;

// Setup callback function for outputs.
    r->heartbeat            = heartbeat;
    r->hash_ctr             = 0;

  // Initial conditions
    for (int i = 0; i < 3; i++) {
        struct reb_particle p = {0};
        p.x = ss_pos[i][0];
        p.y = ss_pos[i][1];
        p.z = ss_pos[i][2];
        p.vx = ss_vel[i][0];
        p.vy = ss_vel[i][1];
        p.vz = ss_vel[i][2];
        p.m = ss_mass[i];
        reb_add(r, p);
    }


    E0 = reb_tools_energy(r);
    L0 = reb_tools_angular_momentum(r);
    printf("E0 = %f",E0);

    reb_integrate(r, tmax);

   // Cleanup
//    reb_mpi_finalize(r);
    reb_free_simulation(r);
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
        int snap_n = (r->hash_ctr);
//        char SNAP[50];
        char SNAP_xyz[50];
//        char SNAP_restart[50];

//    if (reb_output_check(r, 1.0)){  
//        reb_output_timing(r, tmax);
//           }

    if (reb_output_check(r, 0.1)){
//##################################################        
//output of energy error
//##################################################        

	FILE* f = fopen("energy.dat","a");
        reb_integrator_synchronize(r);
        double E = reb_tools_energy(r);  
        double relE = fabs((E-E0)/E0);

        struct reb_vec3d L = reb_tools_angular_momentum(r); 
        struct reb_vec3d relL = {fabs((L.x-L0.x)/L0.x),fabs((L.y-L0.y)/L0.y) , fabs((L.z-L0.z)/L0.z)};    

        fprintf(f,"%e\t %e\t %e\t %e\t %e\t %e\t %e \n",r->t, relE, relL.x, relL.y, relL.z, r->dt, E0);
        fclose(f);

//##################################################        
//output of relative distance
//##################################################        
/*
	FILE* frd= fopen("relative_distance.dat","a");
        reb_integrator_synchronize(r);

        struct reb_particle p0 = r->particles[0];
        struct reb_particle p1 = r->particles[1];
        struct reb_particle p2 = r->particles[2];
	double r12,r23,r31;
        fprintf(of,"%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",i,r->particles[i].m,p.x,p.y,p.z,p.vx,p.vy,p.vz);

        fprintf(frd,"%e\t %e\t %e\t %e\t %e\t %e\t %e \n",r->t, relE, relL.x, relL.y, relL.z, r->dt, E0);
        fclose(frd);
*/
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

