#pragma once


class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp) {
        fscanf(fp, "%lf\n", &time);
        fscanf(fp, "%lld\n", &n_body);
        return n_body;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_body);
    }


};

class FPGrav{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    pot;    

//#######################################################
// add
//#######################################################  
    PS::F64vec jrk;
    PS::F64    ptime;//particle time for individual timestep
    PS::F64    pdt;//particle dt for individual timestep  
    

//for hermite
// ---------------------------------------------------------------------------
    PS::F64vec x0;
    PS::F64vec v0;
    PS::F64vec a0;
    PS::F64vec j0;
    PS::F64vec s0;
    PS::F64vec c0;//position,velocity,acceleration,jerk,snap,crackle
    PS::F64vec xp;
    PS::F64vec vp;
    PS::F64vec a1;
    PS::F64vec j1;


//for collision
// ---------------------------------------------------------------------------

//    std::vector<std::vector<PS::S64>> COL_P;
    static std::vector<PS::S64> COL_P;
//    std::vector<PS::S64> COL_P2;

    static int collisions_N;
    static PS::S32 hermite_step;
    //
    //1:predictor 
    //2:a1j & jerk1j
    //3:new a1j & jerk1j
    //4:next step

    static PS::S64 N_active;
    static PS::S64 id_sun;

    static PS::F64 eps;

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const FPGrav & fp){ 
        mass = fp.mass;
        pos  = fp.pos;
        vel  = fp.vel;

    }

    void copyFromForce(const FPGrav & force) {
        acc = force.acc;
        jrk = force.jrk;
        pot = force.pot;
        COL_P = force.COL_P;
        collisions_N = force.collisions_N;
    }

/*
//for hermite
// ---------------------------------------------------------------------------
    void copyFromForce(const FPGrav & force) {
        a1 = force.a1;
        j1 = force.j1;
        pot = force.pot;
    }
*/
    void clear() {
        acc = 0.0;
        jrk = 0.0;
        pot = 0.0;
//        COL_P.erase(COL_P.begin() , COL_P.end() );
        std::vector<PS::S64>().swap(COL_P);
        COL_P.resize(0);
        collisions_N = 0;
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
        }
};


#ifdef ENABLE_PHANTOM_GRAPE_X86


template <class TParticleJ>
void CalcGravity(const FPGrav * iptcl,
                 const PS::S32 ni,
                 const TParticleJ * jptcl,
                 const PS::S32 nj,
                 FPGrav * force) {
    const PS::S32 nipipe = ni;
    const PS::S32 njpipe = nj;
    PS::F64 (*xi)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64 (*ai)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64  *pi     = (PS::F64  *    )malloc(sizeof(PS::F64) * nipipe);
    PS::F64 (*xj)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * njpipe * PS::DIMENSION);
    PS::F64  *mj     = (PS::F64  *    )malloc(sizeof(PS::F64) * njpipe);
    for(PS::S32 i = 0; i < ni; i++) {
        xi[i][0] = iptcl[i].getPos()[0];
        xi[i][1] = iptcl[i].getPos()[1];
        xi[i][2] = iptcl[i].getPos()[2];
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for(PS::S32 j = 0; j < nj; j++) {
        xj[j][0] = jptcl[j].getPos()[0];
        xj[j][1] = jptcl[j].getPos()[1];
        xj[j][2] = jptcl[j].getPos()[2];
        mj[j]    = jptcl[j].getCharge();
        xj[j][0] = jptcl[j].pos[0];
        xj[j][1] = jptcl[j].pos[1];
        xj[j][2] = jptcl[j].pos[2];
        mj[j]    = jptcl[j].mass;
    }
    PS::S32 devid = PS::Comm::getThreadNum();
    g5_set_xmjMC(devid, 0, nj, xj, mj);
    g5_set_nMC(devid, nj);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);
    for(PS::S32 i = 0; i < ni; i++) {
        force[i].acc[0] += ai[i][0];
        force[i].acc[1] += ai[i][1];
        force[i].acc[2] += ai[i][2];
        force[i].pot    -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}

#else

class CalcGravity{
  public:
  void operator () (const FPGrav * ep_i,
                    const PS::S32 n_ip,
                    const FPGrav * ep_j,
                    const PS::S32 n_jp,
                    FPGrav * force) {
//global const
double M_sun = 1.989*pow(10.0,33.0);//g
double M_earth = 5.9724*pow(10.0,27.0);//g
double one_au=1.49597870*pow(10.0,13.0);//cm
double one_year = 365.0*24.0*60.0*60.0;//year-s
double PI = 3.141592653589793;
double _rho_p = 2.0;// [g/cm^3]
double radius_f = 1.0;//fold enlargement(radius enhancement)

//initializing collsion variable
    force[0].COL_P.erase(force[0].COL_P.begin() , force[0].COL_P.end() );
//    std::vector<PS::S64>().swap(force[0].COL_P);
    force[0].COL_P.resize(0);
    force[0].COL_P.reserve(10);
    force[0].collisions_N = 0;
    int icol[500];
    int jcol[500];
    int col_N = 0;
//    fprintf(stdout, "N_active in CalcGravity : %lld\n", force[0].N_active);
//    fprintf(stdout, "n_jp : %d\n", n_ip);
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;

        force[0].collisions_N = col_N ;

    for(PS::S32 i = 0; i < n_ip; ++i){
      if(i<force[0].N_active){

//        force[0].collisions_N = col_N ;

        PS::F64 mi = ep_i[i].getCharge();
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec vi = ep_i[i].vel;
        PS::F64vec ai = 0.0;
        PS::F64vec jri = 0.0;
        PS::F64 poti = 0.0;
        
        force[i].acc = 0.0;
        force[i].jrk = 0.0;
        force[i].pot = 0.0;
        
        for(PS::S32 j = 0; j < n_jp; ++j){
            if(j != i && j < force[0].N_active  ){
//            if(j != i  ){
            PS::F64vec rij    = xi - ep_j[j].getPos();
            PS::F64vec vij    = vi - ep_j[j].vel;
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);
            PS::F64 mj = ep_j[j].getCharge();

            r3_inv  = r_inv * r_inv * r_inv * mj;
            PS::F64    r5_inv = r_inv * r_inv * r_inv * r_inv * r_inv * mj;
            r_inv  *= mj;

            ai     -= r3_inv * rij;
            PS::F64vec    j1     = r3_inv * vij;
            PS::F64vec    j2     = 3.0 * ( vij * rij ) * r5_inv * rij;
            jri    -= j1 - j2;
            poti   -= r_inv;
            
           

if(force[0].hermite_step == 2 )
{//collision
// ---------------------------------------------------------------------------
            const double mic = mi;
            const double mjc = mj;
            double Rij = sqrt(rij * rij);
            const double rad_i = (3.0*mic*M_sun)/(4.0*PI*_rho_p) ; //unit is AU. ;
            double Ri = radius_f * cbrtl(rad_i)/one_au ; //unit is AU. ;
            const double rad_j = (3.0*mjc*M_sun)/(4.0*PI*_rho_p) ; //unit is AU. ;
            double Rj = radius_f * cbrtl(rad_j)/one_au ; //unit is AU. ;
            double RAD_ij = Ri + Rj;
            if(Rij < RAD_ij){
//            force[0].COL_P.push_back( i );
//            force[0].COL_P.push_back( j );
            icol[force[0].collisions_N] = i ;
            jcol[force[0].collisions_N] = j ;
//            icol.push_back( i );
//            icol.push_back( j );
             ++ force[0].collisions_N ;
//            col_N ++ ;
              }
        }  //if(force[0].hermite_step == 2 )
        
       }
    }
        force[i].acc += ai;
        force[i].jrk += jri;
        force[i].pot += poti;
//        force[i].collisions_N = col_N;
    }
//            fprintf(stdout, "force[0].collisions_N: %lld \n", force[0].collisions_N);
}

//        force[0].collisions_N = col_N;
//        fprintf(stdout, "force[0].collisions_N: %lld \n", force[0].collisions_N);
//        fprintf(stdout, "force[0].collisions_N: %lld \n", col_N);

//        force[0].collisions_N =  col_N ;
}


};

#endif
