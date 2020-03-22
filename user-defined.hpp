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

    static PS::S64 collisions_N;
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
//    force[0].COL_P.resize(1,std::vector<PS::S64>(2));
        force[0].COL_P.erase(force[0].COL_P.begin() , force[0].COL_P.end() );
        std::vector<PS::S64>().swap(force[0].COL_P);
    force[0].COL_P.resize(0);
    force[0].collisions_N = 0;
//    co_n  = 0;
//    fprintf(stdout, "N_active in CalcGravity : %lld\n", force[0].N_active);
//    fprintf(stdout, "n_jp : %d\n", n_ip);
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
      if(i<force[0].N_active){
        PS::F64vec xi = ep_i[i].pos;
        PS::F64vec vi = ep_i[i].vel;
        PS::F64vec ai = 0.0;
        PS::F64vec jri = 0.0;
        PS::F64 poti = 0.0;
        
        force[i].acc = 0.0;
        force[i].jrk = 0.0;
        force[i].pot = 0.0;
        
        for(PS::S32 j = 0; j < n_jp; j++){
            if(j != i && j < force[0].N_active  ){
//            if(j != i  ){
            PS::F64vec rij    = xi - ep_j[j].pos;
            PS::F64vec vij    = vi - ep_j[j].vel;
            PS::F64    r3_inv = rij * rij + eps2;
            PS::F64    r_inv  = 1.0/sqrt(r3_inv);

            r3_inv  = r_inv * r_inv * r_inv * ep_j[j].mass;
            PS::F64    r5_inv = r_inv * r_inv * r_inv * r_inv * r_inv * ep_j[j].mass;
            r_inv  *= ep_j[j].mass;

            ai     -= r3_inv * rij;
            PS::F64vec    j1     = r3_inv * vij;
            PS::F64vec    j2     = 3.0 * ( vij * rij ) * r5_inv * rij;
            jri    -= j1 - j2;
            poti   -= r_inv;

if(force[0].hermite_step == 2 )
{
//collision
// ---------------------------------------------------------------------------
//global const                                                                                                                                           
PS::F64 M_sun = 1.989*pow(10.0,33.0);//g
PS::F64 M_earth = 5.9724*pow(10.0,27.0);//g
PS::F64 one_au=1.49597870*pow(10.0,13.0);//cm
PS::F64 one_year = 365.0*24.0*60.0*60.0;//year-s

            PS::F64 Rij = sqrt(rij * rij);
            PS::F64 _rho_p = 2.0;// [g/cm^3]
//            PS::F64 rho_p = _rho_p * 1.49597871e13 * 1.49597871e13 * 1.49597871e13 / 1.9884e33; // [Msun/AU^3]
            PS::F64 radius_f = 5.0;//fold enlargement(radius enhancement)
            PS::F64 rad_i = radius_f * pow((3.0*ep_i[i].mass*M_sun)/(4.0*M_PI*_rho_p) ,1.0/3.0)/one_au ; //unit is AU. ;
            PS::F64 rad_j = pow((3.0*ep_j[j].mass*M_sun)/(4.0*M_PI*_rho_p) ,1.0/3.0)/one_au * radius_f; //unit is AU. ;
//            PS::F64 RAD_ij = pow(rad_i + rad_j, 2.0);
            PS::F64 RAD_ij = rad_i + rad_j;
            PS::S64 YES = 0;
            
            if(Rij < RAD_ij){
              
            YES = 1 ;

//            fprintf(stdout, "YES: %lld \n", YES);
            for(PS::S64 iC = 0; iC < force[0].collisions_N ; iC++){
            if(force[0].COL_P[2*iC] == j  && force[0].COL_P[2*iC+1] == i || force[0].COL_P[2*iC] == i  && force[0].COL_P[2*iC+1] == j  ){
            YES = 0;
            }
            }

//            fprintf(stdout, "YES: %lld \n", YES);
            if(YES == 1){
  
            force[0].COL_P.push_back( i );
            force[0].COL_P.push_back( j );
//            force[0].COL_P[force[0].collisions_N][0] = i ;
//            force[0].COL_P[force[0].collisions_N][1] = j ;

            fprintf(stdout, "YES: %lld, A collsion occured: %d, %d  \n", YES, i, j);
            force[0].collisions_N ++ ;
//            collisions_N ++ ;
//            col_n = force[0].collisions_N ;
            }

              }

        }
       }
    }

        force[i].acc += ai;
        force[i].jrk += jri;
        force[i].pot += poti;
    }
//            fprintf(stdout, "force[0].collisions_N: %lld \n", force[0].collisions_N);
}
}


};

#endif
