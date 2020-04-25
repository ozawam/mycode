#include<iostream>
#include<fstream>
#include<iomanip>
#include <stdlib.h>
#include<string.h> 
#include <string>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <list>
#include <numeric>
#include <sstream>

using namespace std;

//--------------------
//bar graph collision-n
//--------------------
void c_bar_graph(int bar_end, int n_bin, double c_in, double c_ou, vector<vector<vector<double>>>& c_number, double c_au, int i, int i_a_m, int j_a_m, int Z)
{

   if(c_in <= c_au && c_au < c_ou){
           c_number.at(i_a_m).at(j_a_m).at(Z) +=  1;
// cout <<"a_m.at(i_a_m).at(j_a_m)=" << a_m.at(i_a_m).at(j_a_m) << '\n';
                    };
}


//--------------------
//bar graph e-n
//--------------------
void e_bar_graph(int bar_end, int n_bin, double e_in, double e_ou, vector<vector<double>>& e_number, vector<double> e, int i, int i_a_m, int j_a_m)
{
   if(e_in <= e[i] && e[i] < e_ou){
           e_number.at(i_a_m).at(j_a_m) +=  1;
// cout <<"a_m.at(i_a_m).at(j_a_m)=" << a_m.at(i_a_m).at(j_a_m) << '\n';
                    };
}


//--------------------
//bar graph a-m
//--------------------
void a_m_bar_graph(int bar_end, int n_bin, double a_in, double a_ou, vector<vector<double>>& a_m, vector<double> a_au, vector<double> m, int i, int i_a_m, int j_a_m)
{
/*
if( a_m == nullptr )
    {
        cout << "null" << endl;
        return true;
    }
*/
// cout <<"a[i]=" << a_au[i] << '\n';
   if(a_in <= a_au[i] && a_au[i] < a_ou){
           a_m.at(i_a_m).at(j_a_m) +=  m[i];
// cout <<"a_m.at(i_a_m).at(j_a_m)=" << a_m.at(i_a_m).at(j_a_m) << '\n';
                    }
// return false;
}

/*void a_m_bar_graph_05(double a_in, double a_ou, double a_m_05[15], vector<double> a_au, vector<double> m, int i, int i_a_m)
  {
     a_in += 0.5;
     a_ou += 0.5;
     if(a_in <= a_au[i] && a_au[i] < a_ou){
               a_m_05[i_a_m] += m[i];
  
                      }

  }
  */

//--------------------
//cumulative distribution
//--------------------
double dice(int n,vector<double> value,vector<double> number,double time ,string file_name_cum)
{
  //  vector<double> p(1+n); //確率
  vector<double> f(1+n); //cumulative distribution

  int i;
  int j;
 
    //ファイルへのデータ書き出し

  ofstream cumulative_mass_1(file_name_cum); 
  if(!cumulative_mass_1)
    {
      cout << "I can't open dat file.\n";
      return 1;
    }

  /*
  for(i=1; i<=n ;i++ ){
    p[i]=number[i]/n;
    //    cout << p[i];
    }
  */
  
  for(i=1; i<=n-2 ;i++ )
    {
    for(j=1;j<=i ;j++)
      {
	//f[i] += p[j];
	f[i] += 1.0;
      }
  cumulative_mass_1<< value[i] << "  " << f[i] << "  " << time << '\n';
    }
    return 0;
       
 }    

//--------------------
//main
//--------------------
int main(int argc, char* argv[])
{

  cout <<"You need to delete .dat before starting this analysis." << '\n';
  //単位系cgs
  //無次元化
  double M_sun=1.989*pow(10.0,33.0);//g
  double M_lunar = 7.34581*pow(10.0,25.0);//[g]
  double M_pluto = 1.30*pow(10.0,25.0);//[g]
  double M_earth = 5.972*pow(10.0,27.0);//[g]
  double one_au=1.49597870*pow(10.0,13.0);//cm
  double one_year = 365.0*24.0*60.0*60.0;//year-s
  
  double L_unit=one_au;
  double M_unit=M_sun;
  double G_unit=6.67408*pow(10.0,-8.0);

  double t_16=sqrt(pow(L_unit,3.0)/(G_unit*M_unit));
  double T_unit=t_16;
  double V_unit=L_unit/T_unit;
  //  cout <<"V_unit=" << V_unit << '\n';
//  cout <<"T_unit=" << T_unit << '\n';


  char filename[50];
  char directoryname[100];
  int i;//ファイル内のループ
  int j;//ファイルのループ
  int k;//ディレクトリのループ

  int K;//使うディレクトリ数
  int K_start;
 
  //  K = 9;
  //  cout << "enter number of directories.\n" ;
  //  cin >> K ;
  K_start = (int)strtol(argv[1], NULL, 10);
  K = (int)strtol(argv[2], NULL, 10);

  vector<double> first_line(11);//1行目を読み込ませるための変数
  double time = 0.0 ;
  int time_file_name = 0 ;
  int n; //粒子数

  vector<double> e(100000);//離心率
  vector<double> id(100000);
  vector<double> m(100000);//単位系を直した質量 g
  
  vector<double> mu(100000);//G(m1+m2)
  vector<double> a_au(100000);//軌道長半径、単位AU
  vector<double> a(100000);
  vector<double> omega(100000);//昇交点経度
  vector<double> I(100000);//軌道傾斜角
  vector<double> I_abs(100000);//軌道傾斜角の絶対値
  vector<double> w(100000);//近点引数
  vector<double> p(100000);//半直弦
  vector<double> lambda(100000);
  vector<double> period(100000);
  vector<double> f_ta(100000);
  vector<double> radius(100000);

  vector<double> H(100000);//Hill radius

  vector<double> value(100000); //出る値
  vector<double> number(100000);//その値が出た回数
      
  //相対速度を出すために必要な要素
      
  double v_m;//velocity dispersion m
  vector<double> e_m(K+1);//RMS eccentricity m
  vector<double> I_m(K+1);//RMS inclination m
  
  double r_h ;//hill radius
  double w_kep;//kepler angular velocity
  double v_rel;//relative velocity
  double v_M;//velocity dispersion M
  double v_esc;//escape velocity
  double v_kep;//kepler velocity
  double r_M;//radius M

  double rho = 2.0;//density

  double left,center,right;
  double max,second,min,mean_mass,min_mass,max_mass;
  double asecond,esecond; 
  double sigma_e = 0.0;
  double sigma_I = 0.0;
  double sigma_e_mean = 0.0;
  double sigma_I_mean = 0.0;
  double e_mean;
  double I_mean;
  int n_field;
      
  double f;

  //各効果によるde^2/dt
  double de2_cd_2;
  double de2_cd_1;
  double de2_gas;
  
  double de2_vs_2;
  double de2_vs_1;
  double de2_df;

  double e_2;
  double tcoll22;
  double tcoll21;
  double d2;
  double d1;
  double n2;
  double sigma_coll22;
  double sigma_coll21;
  double rho_pls = 2.0;
  double Sigma2;
  double vesc22;
  double vesc21;
  double F_drag;
  double delta_u;
  double eta = 0.0;
  double rho_gas;
  double Cd = 1.0;
  double Sigma_gas;
  double f_gas ;
   //cout <<"f_gas = ";cin >> f_gas;
  f_gas =  strtod(argv[3], NULL);
  double h_gas ;
  double tau_gas = pow(10.0,6.0);
  double tg22;
  double tg21;
  double sigma_g_22;
  double sigma_g_21;
  double lnA;
  double n1;
  double a_max;
  double e1;
  double a_start_au;//semi major axis[AU]
  double a_start;//semi major axis[cm]

  //  cout <<"a(AU) = ";cin >> a_start_au;
  a_start_au = strtod(argv[4], NULL);
  a_start = a_start_au*one_au;
  
  double da;//width of feeding zone[au]
  double M_target = M_earth; 
  da = 10.0*pow(2.0*M_target/(3.0*M_sun),1.0/3.0)*a_start;
  double width = da/2.0;
  double S2 = 4.0*M_PI*a_start*width;
  double f_ice;
  if(a_start_au <2.7)
    {
      f_ice = 1.0;
    }
  else if(a_start_au >= 2.7)
    {
      f_ice = 4.2;
    }
  double Sigma_solid = f_ice*10.0*pow(a_start_au,-3.0/2.0);

  double k_hcr = 1.41;//heat capacity ratio
  double R_gc = 8.314*100.0*100.0*1000.0 ;//gas constant cm2 g/s2 K mol
  double T = 280.0*pow(a_start_au,-1.0/2.0);//gas temperature
  double M_H2 = 2.0 ;//molecular weight of hydrogen
  //  double c_s = sqrt(k_hcr*R_gc*T/M_H2);//sound velocity cm/s
  double c_s =1.0*pow(10.0,5.0)*sqrt(T/300.0) ;
  
  double M_a; //Mach number
  double c_t = sqrt(8.0/M_PI)*c_s; //the heat velocity
  double R_e; //Reynolds number
  double L_g; //the mean free pass
  double k_B = 1.380658*pow(10.0,-16.0)  ; // Boltzmann constant
  double P_gas; //the Pressure of gas
  double u; //the relative velocity between the gas and particles
  double CD_estimate;
  double ww =0.2;

	  //for bar graph const

	 int i_a_m;
	 int j_a_m;
	 double a_inner;
	 double a_outer;
	 double e_inner;
	 double e_outer;
	 double c_inner;
	 double c_outer;

	 int n_bin = 2;
	 double w_bin = 1.0/n_bin;
	 int bar_start = 0;
	 int bar_end = 30;
	 int e_start = 0;
	 int e_end = 7;
	 int c_start = 1;
	 int c_end = 30;
   int ZMAX = 8+1;
   int Z;
	 vector<vector<double>> a_m(bar_end+1, vector<double>(n_bin)) ;
	 vector<vector<double>> e_number(e_end+1, vector<double>(n_bin)) ;
	 vector<vector<vector<double>>> c_number(c_end+1,vector<vector<double>>(n_bin, vector<double>(ZMAX))) ;

  double c_au;


 //resonance angle
double resonance_angle;
double j_resonance;

 //I/O test
 string buf;

  //データを書き出すファイルのオープン
  ofstream eout("e.dat", ios::app); 
  if(!eout){
    cout << "I can't open e.dat.\n";
    return 1;
  }
      
  ofstream fout("f.dat", ios::app); 
  if(!fout){
    cout << "I can't open f.dat.\n";
    return 1;
  }

  ofstream e_analytical_out("e_analytical.dat", ios::app); 
  if(!e_analytical_out){
    cout << "I can't open e_analytical.dat.\n";
    return 1;
  }
  
  ofstream mout("m.dat", ios::app); 
  if(!mout){
    cout << "I can't open m.dat.\n";
    return 1;
  }
  
  ofstream mout_sigma("m_sigma.dat", ios::app); 
  if(!mout_sigma){
    cout << "I can't open m_sigma.dat.\n";
    return 1;
  }
 
ofstream coout("collision.dat", ios::app); 
   if(!coout){
     cout << "I can't open collision.dat.\n";
     return 1;
   }

/*
ofstream cobout("collision_au_e_before.dat", ios::app); 
   if(!cobout){
     cout << "I can't open collision_au_e_before.dat.\n";
     return 1;
   }
*/

ofstream coaout("collision_au_e_after.dat", ios::app); 
   if(!coaout){
     cout << "I can't open collision_au_e_after.dat.\n";
     return 1;
   }

//--------------------
//CPU TIME modify
//--------------------

ofstream cpuout("cputime.dat", ios::app); 
   if(!coaout){
     cout << "I can't open cputime.dat.\n";
     return 1;
   }
double cpu1,cpu2;
ifstream cpuin("../../cputime.dat" ,ios::in | ios::binary);
while(cpuin >> cpu1 >> cpu2){

cpuout << cpu1 << " " << cpu2  <<'\n';


}

  //collision_after
  string s1 = "merge";
  string s2 = "rebound";
  string s3 = "rebound with crater formation";
  string s4 = "fragmentation";
  
  string filename_co = "../../collision_after.dat";
  string str;
  int countup[5] = {0,0,0,0,0} ;

 ifstream coin(filename_co ,ios::in | ios::binary);
  if(!coin)
  {
       cout << "I can't open the file. \n";
       return 1;
     }

   while(getline(coin,str))
   {
     
  if(equal(str.begin(),str.end(),s1.begin())==true){
    countup[1] ++;
   }

  else if(equal(str.begin(),str.end(),s2.begin())==true){
    countup[2] ++;
    }

 else if(equal(str.begin(),str.end(),s3.begin())==true){
    countup[3] ++;    
     }   

 else if(equal(str.begin(),str.end(),s4.begin())==true){
    countup[4] ++;
     }
        
     }

coout << "merge"  << " " << countup[1] << ","  << "rebound"  << " " << countup[2] << ","  << "rebound with crater formation"  << " " << countup[3] << ","<< "fragmentation"  << " " << countup[4]  <<'\n';
//  fprintf(co,"merge %d\t, rebound %d\t, rebound with crater formation %d\t, fragmentation %d\t\n", countup[1], countup[2], countup[3], countup[4] );
    
    //collision_before

  int N_count = countup[1] + countup[2] + countup[3] + countup[4] ;
  int icb;
//  double before[21];
  double after[24];
  string determination;
//  string filename_cob = "../../collision_before.dat";
  string filename_coa = "../../collision_after.dat";

/*
 ifstream coinb(filename_cob ,ios::in | ios::binary);
  if(!coinb)
  {
       cout << "I can't open the file. \n";
       return 1;
     }
*/

 ifstream coina(filename_coa ,ios::in | ios::binary);
  if(!coina)
  {
       cout << "I can't open the file. \n";
       return 1;
     }

//--------------------
//information of collision
//--------------------

  for(icb=1;icb<=N_count;icb++)
  {
//coinb >> before[0] >> before[1] >>  before[2] >> before[3] >> before[4] >> before[5] >> before[6] >> before[7] >> before[8] >> before[9] >> before[10] >> before[11] >>  before[12] >> before[13] >> before[14] >> before[15] >> before[16] >> before[17] >> before[18] >> before[19] >> before[20] ;

//cobout << before[0]  << " " << before[3] << " " << before[4] << " " << before[13] << " " << before[14] <<'\n';


coina >> after[0] >> after[1] >>  after[2] >> after[3] >> after[4] >> after[5] >> after[6] >> after[7] >> after[8] >> after[9] >> after[10] >> after[11] >>  after[12] >> after[13] >> after[14] >> after[15] >> after[16] >> after[17] >> after[18] >> after[19] >> after[20] >> after[21] >>  after[22] >> after[23] >> after[24] >> after[25] >> after[26] >> after[27] >> after[28] >> after[29] >> after[30] >> after[31] >>  after[32] >> after[33] >>  after[34] >> after[35] >> determination;
//0-3 are state of system
//4-14 are particle 1
//15-25 are particle 2
//26-35 are property of merge particle

coaout << after[0]  << " " << after[6] << " " << after[7] << " " << after[17] << " " << after[18] << " " << after[27] << " " << after[28] <<'\n';

c_au = after[28];

if(after[0] <= 6.0E+3){
    Z = 1;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(6.0E+3 < after[0] && after[0] < 2.5E+4){
    Z = 2;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(2.5E+4 < after[0] && after[0] <  2.0E+5){
    Z = 3;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(2.0E+5 < after[0] && after[0] <  5.0E+5){
    Z = 4;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(5.0E+5 < after[0] && after[0] <  1.0E+6){
    Z = 5;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(1.0E+6 < after[0] && after[0] <  3.7E+6){
    Z = 6;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(3.7E+6 < after[0] && after[0] <  4.8E+6){
    Z = 7;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}
else if(4.8E+6 < after[0] && after[0] <  5.5E+6){
    Z = 8;
	  for(i_a_m =c_start ;i_a_m<=c_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             c_inner = (i_a_m + j_a_m * w_bin)  - w_bin/2.0 ;
             c_outer = (i_a_m + j_a_m * w_bin)  + w_bin/2.0  ;
             c_bar_graph(c_end, n_bin, c_inner, c_outer, c_number, c_au, i, i_a_m, j_a_m, Z);             
    }
  }  
}


for(Z=1 ; Z<ZMAX ;Z++){          
          string file_name_c_n;
          stringstream ss_c_n;

          ss_c_n <<  "c_n_" <<  to_string(Z)  << ".dat";
          ss_c_n >> file_name_c_n;

          ofstream c_nout(file_name_c_n);
          if(!c_nout){
           cout << "I can't open output file.\n";
           return 1;
          }

for(i_a_m =c_start ; i_a_m<=c_end ;i_a_m++){          
          for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
   c_nout << double(i_a_m + j_a_m * w_bin) << " " << c_number.at(i_a_m).at(j_a_m).at(Z)    << '\n';
    }
  }
}

  }//end of information of collision
  
          
  //  for(k=1;k<=K;k++)
for(k=K_start;k<=K;k++)  
  {
      sprintf(directoryname,"../../output/snap%02d",k);
//      cout << directoryname << '\n';
	  
      //読み込み回数の設定
      int N;
      
      int count = 0;
      string name;
      string prefix = "snap";
      struct dirent *de;
      DIR *d = opendir(directoryname);
      while ((de = readdir(d)) != NULL)
	{  
	  name = de -> d_name;
	  //       if (name.find("snap") == 0)
	  if (name.size() >= prefix.size() &&
	      name.compare(0, prefix.size(), prefix) == 0)
	    {
	    
	      count++;
	    }   
	}
      closedir(d);
      printf("ファイルが %d 個あります\n", count);
      
      N = count-1;
      //Nは読み込むファイルの数

    for(j=0;j<=N;j++)
	{
		if(j%10 == 0 || j <= 1000  )
		{
	  //snapshot for animation
	  string file_name_pp;
	  stringstream ss_pp;

	  string file_name_ps;
	  stringstream ss_ps;

	  ss_pp <<  to_string(k)  << "snap_pp" << to_string(j) << ".dat";
	  ss_pp >> file_name_pp;
//	  cout << file_name_pp << endl;

	  ss_ps <<  to_string(k)  << "snap_ps" << to_string(j) << ".dat";
	  ss_ps >> file_name_ps;
//	  cout << file_name_ps << endl;

	  ofstream fout_pp(file_name_pp);
	  if(!fout_pp){
	    cout << "I can't open output file.\n";
	    return 1;
	  }


	  ofstream fout_ps(file_name_ps);
	  if(!fout_ps){
	    cout << "I can't open output file.\n";
	    return 1;
	  }

          //for a_m bar graph 
          string file_name_a_m;
          stringstream ss_a_m;

          ss_a_m <<  "a_m_" << to_string(k) << "_" << to_string(j)  << ".dat";
          ss_a_m >> file_name_a_m;
//          cout << file_name_a_m << endl;

          ofstream a_mout(file_name_a_m);
          if(!a_mout){
           cout << "I can't open output file.\n";
           return 1;
          }

          //for e_n bar graph 
          string file_name_e_n;
          stringstream ss_e_n;

          ss_e_n <<  "e_n_" << to_string(k) << "_" << to_string(j)  << ".dat";
          ss_e_n >> file_name_e_n;
//          cout << file_name_a_m << endl;

          ofstream e_nout(file_name_e_n);
          if(!e_nout){
           cout << "I can't open output file.\n";
           return 1;
          }



          // input of file
	  sprintf(filename,"../../output/snap%02d/snap%05d.dat",k,j);//結果のファイルがどこにあるか確認すること
	  //   cout << filename << '\n';
	  
	  ifstream fin(filename ,ios::in | ios::binary);
	  if(!fin)
	    {
	      cout << "I can't open the file. \n";
	      return 1;
	    }  

	  fin >> n >> first_line[0] >>  first_line[1] >> first_line[2] >> first_line[3] >> first_line[4];

//	   FILE *fp;  /*  ０：ファイルポインタを宣言  */

//	    fp = fopen(filename, "r");       /*  １：ファイルを書き込みモードで開く（オープン）  */

//	    fscanf(fp, "%d %lf %lf %lf %lf %lf", &n, &first_line[0] ,  &first_line[1] , &first_line[2] , &first_line[3] , &first_line[4]);           /*  ２：ファイルから整数を読み出す  */

//	    fclose(fp);                      /*  ３：ファイルを閉じる（クローズ）  */


//	  cout <<"n=" << n << '\n';
//	  cout <<"time(year)=" << first_line[0]  << '\n';
//          cout <<"max(g)=" << first_line[1]  << '\n';
//	  cout <<"min(g)=" << first_line[2]  << '\n';
//          cout <<"RMS e=" << first_line[3]  << '\n';
//	  cout <<"primary mass(Msun)=" << first_line[4]  << '\n';

	  e.resize(n);//離心率
	  id.resize(n);
	  m.resize(n-1);//単位系を直した質量 g
    
	  mu.resize(n);//G(m1+m2)
	  a_au.resize(n);//軌道長半径 
	  a.resize(n);
	  omega.resize(n);//昇交点経度
	  I.resize(n);//軌道傾斜角
	  I_abs.resize(n);//軌道傾斜角の絶対値
	  w.resize(n);//近点引数
	  p.resize(n);//半直弦
	  lambda.resize(n);
	  period.resize(n);
	  f_ta.resize(n);
	  radius.resize(n);

	  H.resize(n);//Hill radius

	  value.resize(n); //出る値
	  number.resize(n);//その値が出た回数

	  n_field = 0;
	  sigma_e = 0.0;
	  sigma_I = 0.0;
	  sigma_e_mean = 0.0;
	  sigma_I_mean = 0.0;

double aJ;
double MJ;
double eJ;
  
	  for(i=0;i<n-1 ;i++)//Because n includes Sun.
	    {
	//      fin >>scientific >> id[i] >> m[i] >>  a_au[i] >> e[i] >> I[i] >> omega[i]  >> w[i] >> lambda[i] >> period[i] >> f_ta[i] >>radius[i];
	fin >>scientific >> id[i] >>scientific >> m[i] >>scientific >>  a_au[i] >>scientific >> e[i] >>scientific >> I[i] >>scientific >> omega[i] >>scientific >> w[i] >>scientific >> lambda[i] >>scientific >> period[i] >>scientific >> f_ta[i] >>scientific >>radius[i];

   //searching Jupiter
 if(m[i] < 1.0*pow(10.0,33.0) && m[i] > 1.0*pow(10.0,30.0)){
   aJ = a_au[i];
   MJ = m[i];
   eJ = e[i];
   }

 getline(fin, buf);
    if (!fin) {
        cerr << "読み込みに失敗" << endl;
        exit(1);
    }

/*
    cout << buf << endl;
 

	  cout <<"time" << id[i] << '\n';
          cout <<"m[g]=" << m[i]  << '\n';
          cout <<"a[au]=" << a_au[i]  << '\n';
          cout <<"e=" << e[i]  << '\n';
          cout <<"I=" << I[i]  << '\n';
          cout <<"Omega=" << omega[i]  << '\n';
          cout <<"omega=" << w[i]  << '\n';
          cout <<"lambda=" << lambda[i]  << '\n';
          cout <<"period=" << period[i]  << '\n';
          cout <<"true anomaly=" << f_ta[i]  << '\n';
          cout <<"radius=" << radius[i]  << '\n';
*/
	      a[i] = a_au[i] * one_au ;
	      v_kep = sqrt(G_unit*M_sun/a[i]);  


//--------------------
//for a-m bar graph and e-n bar graph   
//--------------------

if(m[i] < 1.0e+30){
	  for(i_a_m =bar_start ;i_a_m<=bar_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             a_inner = i_a_m + j_a_m * w_bin ;
             a_outer = i_a_m + (j_a_m + 1.0) * w_bin ;
             a_m_bar_graph(bar_end, n_bin, a_inner, a_outer, a_m, a_au, m, i, i_a_m, j_a_m);             
//             cout <<"a_m =" << a_m.at(i_a_m).at(j_a_m) << '\n';
     //   cout <<"a_m[7]=" << a_m[7] << '\n';

    }
 }  
	  for(i_a_m =e_start ;i_a_m<=e_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             e_inner = (i_a_m + j_a_m * w_bin) * 0.01 ;
             e_outer = (i_a_m + (j_a_m + 1.0) * w_bin) * 0.01;
             e_bar_graph(e_end, n_bin, e_inner, e_outer, e_number, e, i, i_a_m, j_a_m);             
    }
 }  
}
	    }
        // time
          time = first_line[0] ;
          time_file_name = (int)time;

        //resonance angle
        j_resonance = 2.0;
        resonance_angle = ((j_resonance+1.0) - lambda[0] - j_resonance * lambda[10] - w[10])/(2.0 * M_PI); //2πで割ったあまりを出力

// for a-m bar graph
for(i_a_m =bar_start ; i_a_m<=bar_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
 a_mout << double(i_a_m + j_a_m * w_bin) << " " << a_m.at(i_a_m).at(j_a_m) << " "  << time  << "  " <<  n  << '\n';
 a_m.at(i_a_m).at(j_a_m) = 0.0;     
// a_mout << i_a_m+0.5 << " " << a_m_05[i_a_m] << " "  << time  << "  " <<  n  << '\n';
  }
}

// for e-n bar graph
int sum_e_number = 0;
for(i_a_m =e_start ; i_a_m<=e_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
 e_nout << double((i_a_m + j_a_m * w_bin)*0.01) << " " << e_number.at(i_a_m).at(j_a_m) << " "  << time  << "  " <<  n  << '\n';
 sum_e_number += e_number.at(i_a_m).at(j_a_m);
 e_number.at(i_a_m).at(j_a_m) = 0.0;     
// a_mout << i_a_m+0.5 << " " << a_m_05[i_a_m] << " "  << time  << "  " <<  n  << '\n';
  }
}
int e_truncated = n - 2 - sum_e_number  ;//Two is Sun and Jupiter.
e_nout << (e_end + 1.0) * 0.01  << " " << e_truncated  << " "  << time  << "  " <<  n  << '\n';

// cout <<" a_m[7]  = " << a_m[7]  << '\n';
	  //inclinationは右向きに傾いているときに負の値をとるので、絶対値をとる
	  for(i=0 ;i<=n-1 ;i++)
	    {
	      I_abs[i] = abs(I[i]);
	    }
	

	  //output_protoplanet,planetesimals
	  
	  long double pp_max_mass;
	  pp_max_mass = *max_element(m.begin(), m.end());
//	  long double max;
	  max = pp_max_mass;


	if(j == 0){
//	  long double min_mass;
	  min_mass = *min_element(m.begin(), m.end());
	}
//	  long double min;
	  min = min_mass;

	  vector<long double> H(1+n);
	  int Npp =0;

	  //protoplanet_selection
	  //runaway body selection
	  //上：最大質量天体の1/5以上の質量を持っていたら原始惑星とする
	  //下：初期質量(つまり最小天体の質量)の10倍以上の大きさになったらrunaway bodyとする
    //Information of Jupiter is in a header fout_ps.
	  
	fout_ps << first_line[0] << " " << eJ << " " << aJ << " " << MJ << '\n' << '\n' << '\n' ;
	for(i=0 ;i<=n ;i++){
	    /*
	      if(m[i] < pp_max_mass/5.0){
	      fout_ps << e[i] << " " << a[i]/one_au << " " << m[i] << '\n';
	      }
	      else if(m[i] >= pp_max_mass/5.0){
	      H[i]=a[i]*pow(m[i]/(3.0*M_sun),1.0/3.0);
	      fout_pp << e[i] << " " << a[i]/one_au << " " << m[i] << " " << H[i]/one_au << " " << 5.0*H[i]/one_au << '\n';
	      }
	    */

	    if(m[i] < min_mass*10.0){
	      fout_ps << e[i] << " " << a_au[i] << " " << m[i]  << '\n';
	    }
	    else if(m[i] >= min_mass*10.0 && m[i] != MJ){
	      H[i]=a[i]*pow(m[i]/(3.0*M_sun),1.0/3.0);
	      //             fout_pp << e[i] << " " << a[i]/one_au << " " << m[i] << " " << H[i]/one_au << " " << 5.0*H[i]/one_au << '\n';
	      //4行目はポイントサイズ
	      fout_pp << e[i] << " " << a_au[i] << " " << m[i] << " " << m[i]/(pow(10.0,24.0)) << " " <<  first_line[0] << '\n';
	      Npp ++;
	 }
	    
	    
	  }
	  
//--------------------
//second largest mass
//--------------------

	second = 0.0;

	  for(i=0; i<=n-1; i++)
	    {
	      if(second < m[i] && m[i] != max)
		{
		  second = m[i];
		}
}
//--------------------
//Property of max mass
//--------------------

	  for(i=0; i<=n-1; i++)
	    {
	    if(m[i] == second)
		{
		  asecond = a_au[i];
		  esecond = e[i];
		}

	    }

	  
	  mean_mass =( accumulate(&m[0], &m[n-1],0.0) - MJ )/n  ;//mean_mass exclude Jupiter and Sun
//	  cout <<" mean_mass  = " << mean_mass  << '\n'; 
	  mout << first_line[0]  << "  " << max << "  " << mean_mass  << "  " << second  << "  " << asecond << "  " <<  esecond  << "  " << resonance_angle  << "  " <<  n   << "  " << Npp  <<'\n';//g
	  mout_sigma << first_line[0]/pow(a_start_au,3.0/2.0)*Sigma_solid << "  " << max << "  " << mean_mass  << "  " <<  n  << '\n';//g //time[Tkep sigma_solid]

	  // cout <<"min=" << min << '\n';
	  // cout <<"n=" << n << '\n';




	  //sort and cumulative mass
	                        
   string file_name_cum;
   stringstream ss_cum;
                     
  //▸-  ss_cum <<  to_string(k)  << "cumulative" << to_string(j) << ".dat";
      ss_cum <<  "cumulative" << to_string(j) << ".dat";
   ss_cum >> file_name_cum;
//   cout << file_name_cum << endl;
                    
   for(i=0 ;i<=n-2 ;i++)
     {            
       value[i] = m[i] ;
       number[i] = 1.0;
     }            
   stable_sort(value.begin(),value.end(),greater<double>());
 //▸-  time = first_line[0]*(T_unit/one_year) ;
 //        time_file_name = (int)time;
   dice(n,value,number,time,file_name_cum);
	      
	  for(i=0;i<=n-1;i++)
	    {
	      if( k == 1 && j == 0 )
		{
		  n_field = n;
		  v_kep = sqrt(G_unit*M_sun/a_start);
		  v_M = sqrt(pow(e[i],2.0)+pow(I[i],2.0))*v_kep;
		  r_M = pow(3.0*m[i]/(4.0*M_PI*rho),1.0/3.0) ;
		  v_esc = sqrt(2.0*G_unit*m[i]/r_M);
		  r_h = pow(2.0*m[i]/(3.0*M_sun),1.0/3.0)*a[i];
		  a_max = a_start;
		  e1 = e[i];
		  //	  w_kep = sqrt(G_unit*M_sun/pow(a[i],3.0));
		  left = 2.0*r_h*w_kep;
		  right = v_esc;
		  sigma_e += pow(e[i],2.0);
		  sigma_I += pow(I[i],2.0); 
 
		  sigma_e_mean += e[i];
		  sigma_I_mean += I[i]; 
		}
	      else
		{
// Sigma of e and i
	if(m[i] < 1.0e+30){
		sigma_e += pow(e[i],2.0);
		sigma_I += pow(I[i],2.0); 
	
		sigma_e_mean += e[i];
		sigma_I_mean += I[i]; 
		      
		n_field = n;
			}
      
		  if(m[i] == max)
		    {
		      v_kep = sqrt(G_unit*M_sun/a[i]);
		      v_M = sqrt(pow(e[i],2.0)+pow(I[i],2.0))*v_kep;
		      r_M = pow(3.0*m[i]/(4.0*M_PI*rho),1.0/3.0) ;
		      v_esc = sqrt(2.0*G_unit*m[i]/r_M);
		      r_h = pow(2.0*m[i]/(3.0*M_sun),1.0/3.0)*a[i];
		      a_max = a[i];
		      e1 = e[i];
		      //	  w_kep = sqrt(G_unit*M_sun/pow(a[i],3.0));
		      left = 2.0*r_h*w_kep;
		      right = v_esc;
		      
		    }
		  else if(m[i] == min)
		    {
//		      sigma_e += pow(e[i],2.0);
//		      sigma_I += pow(I[i],2.0); 
	
//		      n_field += 1;
		      
		    }
		}
	    }
	  /*
	  //all
	   for(i=0;i<=n-1;i++)
	    {
	     
		  sigma_e += pow(e[i],2.0);
		  sigma_I += pow(I[i],2.0); 
		  if(m[i] == max)
		    {
		      v_kep = sqrt(G_unit*M_sun/a[i]);  
		      v_M = sqrt(pow(e[i],2.0)+pow(I[i],2.0))*v_kep;
		      r_M = pow(3.0*m[i]/(4.0*M_PI*rho),1.0/3.0) ;
		      v_esc = sqrt(2.0*G_unit*m[i]/r_M);
		      r_h = pow(m[i]/(3.0*M_sun),1.0/3.0)*a[i];
		      w_k = sqrt(G_unit*M_sun/pow(a[i],3.0));
		      left = 2.0*r_h*w_k;
		      right = v_esc;	  
		    }

	    }
	  */
	  v_kep = sqrt(G_unit*M_sun/a_start);
	  w_kep = sqrt(G_unit*M_sun/pow(a_start,3.0));
	  //	  e_m[k] = sqrt(sigma_e/n);
	  //	  I_m[k] = sqrt(sigma_I/n);
	  e_m[k] = sqrt(sigma_e/n_field);
	  I_m[k] = sqrt(sigma_I/n_field);
	  e_mean = sigma_e/n_field;
	  I_mean = sigma_I/n_field;

	  v_m = sqrt(pow(e_m[k],2.0)+pow(I_m[k],2.0))*v_kep;//random velocity of pls

	  v_rel = sqrt(pow(v_M,2.0)+pow(v_m,2.0));
  
	  center = v_rel;
	  f = 1+ pow(v_esc,2.0)/pow(v_rel,2.0);


	  
	  //analytical solution
	  //collisional damping 2-2
	  d2 = pow((3.0*min)/(4.0*M_PI*rho_pls),1.0/3.0)  ;
	  vesc22 = sqrt(2.0*G_unit*(min+min)/(d2+d2))  ;
	  Sigma2 = min * n_field/S2  ;
	  n2 = Sigma2*w_kep/(min*v_m);
	  sigma_coll22 = M_PI*pow(2.0*d2,2.0)*(1.0+(pow(vesc22,2.0)/pow(v_m,2.0))) ;
	  //   tcoll22 = min /( Sigma2 *w_kep )*M_PI*4.0*pow(d2,2.0)*(1+(pow(vesc22,2.0))/(pow(v_m,2.0))) ;
	  tcoll22 = 1.0/(n2*sigma_coll22*v_m)  ;
	  de2_cd_2 = 1.0/2.0 * pow(e_m[k],2.0)/tcoll22  ;

	  //collisional damping 2-1
	  d1 = pow((3.0*max)/(4.0*M_PI*rho_pls),1.0/3.0)  ;
	  vesc21 = sqrt(2.0*G_unit*(max+min)/(d1+d2))  ;
	  n1 = 1.0/(2.0*M_PI*a_max*10.0*r_h)*w_kep/(v_m);
	  sigma_coll21 = M_PI*pow(d1+d2,2.0)*(1.0+(pow(vesc21,2.0)/pow(v_m,2.0))) ;
	  //   tcoll22 = min /( Sigma2 *w_kep )*M_PI*4.0*pow(d2,2.0)*(1+(pow(vesc22,2.0))/(pow(v_m,2.0))) ;
	  tcoll21 = 1.0/(n1*sigma_coll21*v_m)  ;
	  de2_cd_1 = pow(e_m[k],2.0)/tcoll21  ;	  

	  //gas drag
	  delta_u = sqrt(5.0/8.0 * pow(e_m[k],2.0) + 1.0/2.0 *  pow(1.0/2.0*e_m[k],2.0) +eta )* v_kep;
	  Sigma_gas = 2400.0 * f_gas *pow(a_start_au,-3.0/2.0)*exp(-1.0*first_line[0]/tau_gas) ;
	  h_gas = 0.047*pow(a_start_au,5.0/4.0)*one_au;
	  rho_gas = Sigma_gas/h_gas ;
	  F_drag = 1.0/2.0 * Cd *M_PI * pow(d2,2.0) * rho_gas *pow(delta_u,2.0)  ; 
	  de2_gas = 2.0* pow(e_m[k],2.0) /(min * delta_u) *F_drag  ;

	  //viscous stirring 2-2
	  lnA = 5.0;
	  sigma_g_22 =  M_PI * pow((G_unit*(2.0*min)/(pow(v_m,2.0))),2.0)*lnA;
	  tg22 = 1.0/(n2 * sigma_g_22 * v_m) ;
	  de2_vs_2 =  pow(e_m[k],2.0)/(2.0*tg22) ;

	  //viscous stirring 2-1
	  sigma_g_21 =  M_PI * pow((G_unit*(min+max)/(pow(v_m,2.0))),2.0)*lnA;
	  //	  n1 = 1.0/(2.0*M_PI*a_max*10.0*r_h)*w_kep/(v_m);
	  tg21 = 1.0/(n1 * sigma_g_21 * v_m) ;
	  de2_vs_1 =  pow(e_m[k],2.0)/(tg21) ;

	  //dynamical friction
	  de2_df =  max/(max+min) * (max*pow(e1,2.0) - min*pow(e_m[k],2.0))/(max+min) /tg21;

	  //Mach number
	  u = sqrt(3.0)/2.0 *e_m[k] *v_kep;
	  M_a = u/c_t;

	  //Reynolds number
	  P_gas = pow(c_s,2.0) *rho_gas  ;
	  L_g = 1.0/(sqrt(2.0)*M_PI*pow(2.0*d2,2.0)) * k_B * T/P_gas ;
	  R_e = (6.0*d2*u)/(L_g*c_t);

	  //estimate of CD
	  CD_estimate = ((2.0-ww)*M_a)/(1.0+M_a) + ww  ;
	  
	  //  cout <<"left = " << left  << '\n';
	  //  cout <<"center = " << center  << '\n';
	  //  cout <<"right = " << right  << '\n';
	  //  cout <<"unit cm/s " << '\n';
	  /*
	  cout <<" d2 = " << d2  << '\n';
	  cout <<" vesc22 = " << vesc22  << '\n';
	  cout <<" Sigma2 = " << Sigma2  << '\n';
	  cout <<" n_field = " << n_field  << '\n';
	  */

	  /*
	  cout <<" The sound velocity  = " << c_s  << '\n';
	  cout <<" Mach number = " << M_a  << '\n';
	  cout <<" Reynolds number  = " << R_e  << '\n';
	  cout <<" CD estimate  = " << CD_estimate  << '\n';
	  
	 
	  cout <<" de2_cd_2 = " << de2_cd_2 << '\n';
	  cout <<" de2_cd_1 = " << de2_cd_1 << '\n';
	  cout <<" de2_gas  = " << de2_gas  << '\n';
	  cout <<" de2_vs_2  = " << de2_vs_2  << '\n';
	  cout <<" de2_vs_1  = " << de2_vs_1  << '\n';
	  cout <<" de2_df  = " << de2_df  << '\n';
          */
    
	  eout << first_line[0] << "  " <<  e_m[k] << "  " << e_mean  <<  '\n';
	  fout << first_line[0] << "  " << f << '\n';
	  e_analytical_out << first_line[0] << "  " << de2_cd_1 << "  "<< de2_cd_2 << "  " << de2_gas << "  " << de2_vs_1 << "  " << de2_vs_2 << "  " << de2_df <<  '\n';
          cout << resetiosflags(ios_base::floatfield);
  //     fclose(fp);                      /*  ３：ファイルを閉じる（クローズ）  */
	}

	}

    }	
 
      return 0;

    
}


