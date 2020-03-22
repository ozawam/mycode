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
#include "integrator.hpp"

PS::F64 FPGrav::eps = 0.00;//1.0/64.0 ;//1.0e-4;
std::vector<PS::S64> FPGrav::COL_P;
PS::S64 FPGrav::collisions_N;
PS::S32 FPGrav::hermite_step;
PS::S64 FPGrav::N_active;
PS::S64 FPGrav::id_sun;

template<class Tpsys>
void calcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
    if(clear){
        etot = ekin = epot = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
//        epot_loc += system[i].mass * (system[i].pot + system[i].mass / FPGrav::eps);
       epot_loc += system[i].mass * (system[i].pot );
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    etot_loc  = ekin_loc + epot_loc;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);    
}

void printHelp() {
    std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"s: time_step (default: 1.0 / 128.0)"<<std::endl;
    std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
    std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    PS::S32 ret;
    char dir1[1024];
    char dir2[1024];
    
    sprintf(dir1, "%s/original", dir_name);
    sprintf(dir2, "%s/snap", dir_name);

    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
            mkdir("./result/original/",0777);//include N_active
            mkdir("./reslut/snap/",0777);//only N_active
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}


int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

//#######################################################
// cpu time
//#######################################################  
   struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
//  printf("time:    %10ld.%09ld CLOCK_REALTIME\n", ts.tv_sec, ts.tv_nsec); 
  double time_before = ts.tv_sec+(ts.tv_nsec/1000000000.0) ;

    PS::Initialize(argc, argv);
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 time_end = 10.0;
    PS::F64 dt = 1.0 / 128.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr << "theta =" << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step = " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit = " << n_group_limit << std::endl;
            break;
        case 'N':
            n_tot = atoi(optarg);
            std::cerr << "n_tot = " << n_tot << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
     fprintf(stdout, "CHECK0 \n");
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }

     fprintf(stdout, "CHECK1 \n");

//#######################################################
// set of output and time end
//#######################################################  

   	//constant
        PS::F64 M_sun=1.989*pow(10.0,33.0);//g
        PS::F64 one_au=1.49597870*pow(10.0,13.0);//cm
        PS::F64 one_year = 365.0*24.0*60.0*60.0;//year-s    
        PS::F64 L_unit=one_au;
        PS::F64 M_unit=M_sun;
        PS::F64 G_unit= 6.67408*pow(10.0,-8.0);

        PS::F64 t_16=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));
        PS::F64 T_unit=t_16;

//conversion of time unit. year->year
        dt_snap = dt_snap/(T_unit/one_year);
        time_end = time_end/(T_unit/one_year) ;
//        time_end = time_end/0.15924 ;

        fprintf(stdout, "T_unit: %e time_end: %e  \n",
                        T_unit/one_year, time_end);
     fprintf(stdout, "CHECK2 \n");
    makeOutputDirectory(dir_name);

    std::ofstream fout_eng;

    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
	sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    std::ofstream fout_cpu;

    if(PS::Comm::getRank() == 0) {
        char sout_cpu[1024];
	sprintf(sout_cpu, "%s/cpu_time.dat", dir_name);
        fout_cpu.open(sout_cpu);
    }
//#######################################################
// initial condition
//#######################################################  

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    PS::S32 n_loc    = 0;
    PS::F32 time_sys = 0.0;
    if(PS::Comm::getRank() == 0) {
	FileHeader header;
//	system_grav.readParticleAscii("INITIAL_CONDITION/output/snap00000.dat", header);
//	system_grav.readParticleAscii("INITIAL_CONDITION_ONE_KEPLER/output/snap00000.dat", header);
//	system_grav.readParticleAscii("INITIAL_CONDITION_PYTHAGORAS/output/snap00000.dat", header);
//	system_grav.readParticleAscii("INITIAL_CONDITION_KOKUBO_IDA_1996/output/snap00000.dat", header);
//	system_grav.readParticleAscii("INITIAL_CONDITION_KOKUBO_IDA_1996N300/output/snap00000.dat", header);
//	system_grav.readParticleAscii("INITIAL_CONDITION_KOKUBO_IDA_1996N50/output/snap00000.dat", header);
	system_grav.readParticleAscii("INITIAL_CONDITION_BEAUGE_ten_particles//output/snap00000.dat", header);
	time_sys = header.time;
    system_grav[0].N_active =  system_grav.getNumberOfParticleGlobal();
    fprintf(stdout, "N_active: %lld\n", system_grav[0].N_active);

    } else {
        system_grav.setNumberOfParticleLocal(n_loc);
    }

    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_grav);
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif
    PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

    const PS::S32 n_walk_limit = 200;
    const PS::S32 tag_max = 1;
#ifdef MULTI_WALK
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                RetrieveKernel,
                                                tag_max,
                                                system_grav,
                                                dinfo,
                                                n_walk_limit);
#else
    tree_grav.calcForceAllAndWriteBack(CalcGravity(),
                                       system_grav,
                                       dinfo);
#endif
    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    calcEnergy(system_grav, Etot0, Ekin0, Epot0);
    
    fprintf(stdout, "initial energy: %+e\n", Etot0);

    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;
    PS::TimeProfile TP ;
    PS::F64 cpu[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0} ;
    PS::F64 cputime_output = 0.0;
    PS::F64 time_output_last = 0.0;
    PS::F64 time_output_now = 0.0;

    initial_timestep( system_grav, dt, time_sys);

//#######################################################
// start of calculation
//#######################################################  

    while(time_sys < time_end){

//#######################################################
// output snap and cputime
//#######################################################  
        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
            char filename[256];
            sprintf(filename, "%s/original/snap%05d.dat", dir_name, id_snap);
            FileHeader header;
            header.time   = time_sys;
            header.n_body = system_grav.getNumberOfParticleGlobal();
            system_grav.writeParticleAscii(filename, header);
            time_snap += dt_snap;
        
            char filename_snap[256];
            sprintf(filename_snap, "%s/snap/snap%05d.dat", dir_name, id_snap++);
            std::ofstream fout;
            fout.open(filename_snap,std::ios::out);
    
            fout << time_sys << std::endl;
            fout << system_grav[0].N_active << std::endl;
            
            PS::F64vec pos_s =system_grav[system_grav[0].id_sun].pos;
            PS::F64vec vel_s =system_grav[system_grav[0].id_sun].vel;
            PS::F64vec pos_c_sun;
            PS::F64vec vel_c_sun;
            
            for(PS::S64 J=0; J<system_grav[0].N_active; J++){
            pos_c_sun = system_grav[J].pos - pos_s;
            vel_c_sun = system_grav[J].vel - vel_s;

            fout << J << "  " << system_grav[J].mass
                      << "  " << pos_c_sun[0]  <<  "  " <<  pos_c_sun[1]  << "  " <<  pos_c_sun[2]
                     << "  " << vel_c_sun[0]  <<  "  " <<  vel_c_sun[1]  << "  " << vel_c_sun[2]
                     <<  std::endl;
            }
            
            calcEnergy(system_grav, Etot1, Ekin1, Epot1);
        if(PS::Comm::getRank() == 0){
            TP = dinfo.PS::DomainInfo::getTimeProfile();
            TP = system_grav.getTimeProfile();
            TP = tree_grav.getTimeProfile();

            fout_eng << time_sys << "  " << fabs((Etot1 - Etot0) / Etot0)  << "  " << dt  <<  std::endl;
            //unit of cpu time is second
            clock_gettime(CLOCK_REALTIME, &ts);
            time_output_now = ts.tv_sec+(ts.tv_nsec/1000000000.0) ;
            cputime_output  = time_output_now - time_output_last ;
            
            fout_cpu << time_sys 
                     << "  " << cputime_output
                     << "  " << cpu[0]
                     << "  " << cpu[1]  <<  "  " << cpu[2] 
                     << "  " << cpu[3] 
                     << "  " << cpu[4]  <<  "  " << cpu[5]
                     << "  " << cpu[6]  <<  "  " << cpu[7]
                     << "  " << cpu[8]  <<  "  " << cpu[9]
                     << "  " << cpu[10]  <<  "  " << cpu[11]  <<  "  " << cpu[12]
                     << "  " <<  std::endl;

            cpu[0] = 0.0 ;
            cpu[1] = 0.0 ;
            cpu[2] = 0.0 ;
            cpu[3] = 0.0 ;
            cpu[4] = 0.0 ;
            cpu[5] = 0.0 ;
            cpu[6] = 0.0 ;
            cpu[7] = 0.0 ;
            cpu[8] = 0.0 ;
            cpu[9] = 0.0 ;
            cpu[10] = 0.0 ;
            cpu[11] = 0.0 ;
            cpu[12] = 0.0 ;
            time_output_last = time_output_now ;

            TP.PS::TimeProfile::clear();
            dinfo.clearTimeProfile();
            system_grav.clearTimeProfile();
            tree_grav.clearTimeProfile();
          // In log.log
          //  fprintf(stdout, "time: %10.7f energy error: %+e timestep: %10.7f \n",
          //              time_sys, fabs((Etot1 - Etot0) / Etot0), dt);
        }
        }

        
//#######################################################
// orbit integration
//#######################################################  
/*
leap_frog(system_grav,
          dt,
          time_sys,
          n_loop,
          dinfo,
          tag_max,
          n_walk_limit,
          tree_grav
);
*/

  hermite(system_grav,
          dt,
          time_sys,
          n_loop,
          dinfo,
          tag_max,
          n_walk_limit,
          tree_grav
);
        n_loop++;

        cpu[0] += TP.PS::TimeProfile::getTotalTime() * 1.0e-3  ;
        cpu[1] += TP.collect_sample_particle * 1.0e-3   ;
        cpu[2] += TP.decompose_domain * 1.0e-3    ;
        cpu[3] += TP.exchange_particle * 1.0e-3   ;
        cpu[4] += TP.make_local_tree * 1.0e-3  ;
        cpu[5] += TP.make_global_tree * 1.0e-3   ;
        cpu[6] += TP.calc_force * 1.0e-3  ;
        cpu[7] += TP.calc_moment_local_tree * 1.0e-3  ;
        cpu[8] += TP.calc_moment_global_tree  ;
        cpu[9] += TP.make_LET_1st * 1.0e-3   ;
        cpu[10] += TP.make_LET_2nd * 1.0e-3  ;
        cpu[11] += TP.exchange_LET_1st * 1.0e-3   ;
        cpu[12] += TP.exchange_LET_2nd * 1.0e-3  ;

        TP.PS::TimeProfile::clear();
        dinfo.clearTimeProfile();
        system_grav.clearTimeProfile();
        tree_grav.clearTimeProfile();

    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

//#######################################################
// cpu time
//#######################################################  
  clock_gettime(CLOCK_REALTIME, &ts);
  double time_after = ts.tv_sec+(ts.tv_nsec/1000000000.0) ;
         std::cout << "cputime:" << "  " << time_after - time_before ;
    PS::Finalize();
    return 0;
}
