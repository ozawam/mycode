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

  char filename[50];
  char directoryname[100];
  int i;//ファイル内のループ
  int j;//ファイルのループ
  int k;//ディレクトリのループ

  int K;//使うディレクトリ数
  int K_start;
 
  K_start = (int)strtol(argv[1], NULL, 10);
  K = (int)strtol(argv[2], NULL, 10);

  double time = 0.0 ;
  double time_year;
  int time_file_name = 0 ;
  int n; //粒子数

  vector<double> e(100000);//離心率
  vector<int> id(100000);
  vector<double> m_k(100000);//規格化された質量 1M_sun
  vector<double> m(100000);//単位系を直した質量 g

  //単位系を直した位置 cm
  vector<double> x(100000);
  vector<double> y(100000);
  vector<double> z(100000);
  vector<double> r(100000);//中心からの距離 m
      
  //読み込む位置 au
  vector<double> x_au(100000);
  vector<double> y_au(100000);
  vector<double> z_au(100000);
  
  //単位系を直した速度 cm/s
  vector<double> v_x(100000);
  vector<double> v_y(100000);
  vector<double> v_z(100000);

  //読み込む速度 au/T_unit
  vector<double> v_x_cm(100000);
  vector<double> v_y_cm(100000);
  vector<double> v_z_cm(100000);

  //速度ベクトルの絶対値
  vector<double> V(100000); 
  
  vector<double> mu(100000);//G(m1+m2)
  vector<double> a(100000);//軌道長半径 
  vector<double> h(100000);//角運動量
  vector<double> hx(100000);//角運動量の3成分
  vector<double> hy(100000);//角運動量の3成分
  vector<double> hz(100000);//角運動量の3成分
  vector<double> omega(100000);//昇交点経度
  vector<double> I(100000);//軌道傾斜角
  vector<double> I_abs(100000);//軌道傾斜角の絶対値
  vector<double> w(100000);//近点引数
  vector<double> p(100000);//半直弦

  vector<double> H(100000);//Hill radius
  
  vector<double> value(100000); //出る値
  vector<double> number(100000);//その値が出た回数
      
//  vector<double> e_m(K+1);//RMS eccentricity m
//  vector<double> I_m(K+1);//RMS inclination m
  
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
  double sigma_e_rms = 0.0;
  double sigma_I_rms = 0.0;
  double e_rms;
  double I_rms;
  double sigma_e_mean = 0.0;
  double sigma_I_mean = 0.0;
  double e_mean;
  double I_mean;
  double sigma_e_rms_hill = 0.0;
  double sigma_I_rms_hill = 0.0;
  double e_rms_hill;
  double I_rms_hill;
  int n_field;
      
  double f;

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

/*
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
*/
/*
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

 e
 
 lse if(equal(str.begin(),str.end(),s3.begin())==true){
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


// ifstream coinb(filename_cob ,ios::in | ios::binary);
//  if(!coinb)
//  {
//       cout << "I can't open the file. \n";
//       return 1;
//     }


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
*/  
          
//--------------------
//analysis of snap
//--------------------
for(k=K_start;k<=K;k++)  
  {
      sprintf(directoryname,"../result");
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
		if(j%10 == 0  )
		{
	  //snapshot for animation
//--------------------
//set output
//--------------------
  //データを書き出すファイルのオープン
  ofstream eout("e.dat", ios::app); 
  if(!eout){
    cout << "I can't open e.dat.\n";
    return 1;
  }
      
  
  ofstream mout("m.dat", ios::app); 
  if(!mout){
    cout << "I can't open m.dat.\n";
    return 1;
  }
  
ofstream coout("collision.dat", ios::app); 
   if(!coout){
     cout << "I can't open collision.dat.\n";
     return 1;
   }


ofstream coaout("collision_au_e_after.dat", ios::app); 
   if(!coaout){
     cout << "I can't open collision_au_e_after.dat.\n";
     return 1;
   }

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
//--------------------
//input dat
//--------------------
	  sprintf(filename,"../result/snap%05d.dat",j);//結果のファイルがどこにあるか確認すること
	  //   cout << filename << '\n';
	  
	  ifstream fin(filename ,ios::in | ios::binary);
	  if(!fin)
	    {
	      cout << "I can't open the file. \n";
	      return 1;
	    }  

	  fin >> time  >> n ;
//	  cout <<"n=" << n << '\n';
//	  cout <<"time(year)=" << first_line[0]*(T_unit/one_year)  << '\n';
       
//	  cout <<"n=" << n << '\n';
	  e.resize(n);//離心率
	  id.resize(n);
	  m_k.resize(n);//規格化された質量 1M_sun
	  m.resize(n);//単位系を直した質量 g
    
	  //単位系を直した位置 cm
	  x.resize(n);
	  y.resize(n);
	  z.resize(n);
	  r.resize(n);//中心からの距離 m
    
	  //読み込む位置 au
	  x_au.resize(n);
	  y_au.resize(n);
	  z_au.resize(n);
    
	  //単位系を直した速度 cm/s
	  v_x.resize(n);
	  v_y.resize(n);
	  v_z.resize(n);

	  //読み込む速度 T_unit cm/s
	  v_x_cm.resize(n);
	  v_y_cm.resize(n);
	  v_z_cm.resize(n);
	  
	  //速度ベクトルの絶対値
	  V.resize(n); 
	  
	  mu.resize(n);//G(m1+m2)
	  a.resize(n);//軌道長半径 
	  h.resize(n);//角運動量
	  hx.resize(n);//角運動量の3成分
	  hy.resize(n);//角運動量の3成分
	  hz.resize(n);//角運動量の3成分
	  omega.resize(n);//昇交点経度
	  I.resize(n);//軌道傾斜角
	  I_abs.resize(n);//軌道傾斜角の絶対値
	  w.resize(n);//近点引数
	  p.resize(n);//半直弦
	  
	  H.resize(n);//Hill radius


	  value.resize(n); //出る値
	  number.resize(n);//その値が出た回数

	  n_field = 0;
	  sigma_e_rms = 0.0;
	  sigma_I_rms = 0.0;
	  sigma_e_mean = 0.0;
	  sigma_I_mean = 0.0;
    sigma_e_rms_hill = 0.0;
    sigma_I_rms_hill = 0.0;

double aJ = 0.0;
double MJ = 0.0;
double eJ = 0.0;
  
	  for(i=0;i<=n-1 ;i++)//Because n includes Sun.
	    {
        //Unit is time T_unit~0.15925 year, mass Msun, length AU, G=1.  
	      fin >> id[i] >> m_k[i] >>  x[i] >> y[i] >> z[i] >> v_x[i] >> v_y[i] >> v_z[i] ;
/*	      //単位系の換算
	      x[i]=x_au[i]*L_unit;
	      y[i]=y_au[i]*L_unit;
	      z[i]=z_au[i]*L_unit;
    
	      v_x[i]=v_x_cm[i]*V_unit; //*30*pow(10.0,3.0);
	      v_y[i]=v_y_cm[i]*V_unit; //*30*pow(10.0,3.0);
	      v_z[i]=v_z_cm[i]*V_unit; //*30*pow(10.0,3.0);
*/    
	      m[i]= m_k[i]*M_unit; //unit:g
    
	      //軌道要素の計算
//	      mu[i]=G_unit*(1.0*M_unit+m[i]);//G*(Msun+m)
	      mu[i]=1.0*(1.0+m_k[i]);
  
	      r[i]=sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));
	      V[i]=sqrt(pow(v_x[i],2)+pow(v_y[i],2)+pow(v_z[i],2));
    
	      hx[i]=y[i]*v_z[i]-z[i]*v_y[i];
	      hy[i]=z[i]*v_x[i]-x[i]*v_z[i];
	      hz[i]=x[i]*v_y[i]-y[i]*v_x[i];
	      h[i]=sqrt(pow(hx[i],2.0)+pow(hy[i],2.0)+pow(hz[i],2.0));
	      //  h[i]=sqrt(pow(y[i]*v_z[i]-z[i]*v_y[i],2.0)+pow(z[i]*v_x[i]-x[i]*v_z[i],2.0)+pow(x[i]*v_y[i]-y[i]*v_x[i],2.0));

	      p[i]=pow(h[i],2.0)/mu[i];
	      a[i]=1.0/((2.0/r[i])-(pow(V[i],2.0)/mu[i]));
	    
	      e[i]=sqrt(1.0-(p[i]/a[i]));

	      I[i]=atan(sqrt(pow(hx[i],2.0)+pow(hy[i],2.0))/hz[i]);
	      //    I[i]=atan(sqrt(pow(y[i]*v_z[i]-z[i]*v_y[i],2.0)+pow(z[i]*v_x[i]-x[i]*v_z[i],2.0))/(x[i]*v_y[i]-y[i]*v_x[i]));
//	      v_kep = sqrt(G_unit*M_sun/a[i]);  
	      v_kep = sqrt(1.0*1.0/a[i]);  

	      H[i]=a[i]*pow(m_k[i]/(3.0*1.0),1.0/3.0);
	   
   //searching Jupiter
 if(m[i] < 1.0*pow(10.0,33.0) && m[i] > 1.0*pow(10.0,30.0)){
   aJ = a[i];
   MJ = m[i];
   eJ = e[i];
   }

//--------------------
//for a-m bar graph and e-n bar graph   
//--------------------

if(m[i] < 1.0e+30){
	  for(i_a_m =bar_start ;i_a_m<=bar_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
             a_inner = i_a_m + j_a_m * w_bin ;
             a_outer = i_a_m + (j_a_m + 1.0) * w_bin ;
             a_m_bar_graph(bar_end, n_bin, a_inner, a_outer, a_m, a, m, i, i_a_m, j_a_m);             
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
//        cout <<"H[7]=" << H[7] << '\n';
        // time
          time_year = time * (T_unit/one_year) ;
          time_file_name = (int)time_year;
/*
        //resonance angle
        j_resonance = 2.0;
        resonance_angle = ((j_resonance+1.0) - lambda[0] - j_resonance * lambda[10] - w[10])/(2.0 * M_PI); //2πで割ったあまりを出力
*/
// for a-m bar graph
for(i_a_m =bar_start ; i_a_m<=bar_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
 a_mout << double(i_a_m + j_a_m * w_bin) << " " << a_m.at(i_a_m).at(j_a_m) << " "  << time_year  << "  " <<  n  << '\n';
 a_m.at(i_a_m).at(j_a_m) = 0.0;     
  }
}

// for e-n bar graph
int sum_e_number = 0;
for(i_a_m =e_start ; i_a_m<=e_end ;i_a_m++){
        for(j_a_m =0; j_a_m<n_bin ;j_a_m++){
 e_nout << double((i_a_m + j_a_m * w_bin)*0.01) << " " << e_number.at(i_a_m).at(j_a_m) << " "  << time_year  << "  " <<  n  << '\n';
 sum_e_number += e_number.at(i_a_m).at(j_a_m);
 e_number.at(i_a_m).at(j_a_m) = 0.0;     
  }
}
int e_truncated = n - 2 - sum_e_number  ;//Two is Sun and Jupiter.
e_nout << (e_end + 1.0) * 0.01  << " " << e_truncated  << " "  << time_year  << "  " <<  n  << '\n';

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

	  int Npp =0;

	  //protoplanet_selection
	  //runaway body selection
	  //上：最大質量天体の1/5以上の質量を持っていたら原始惑星とする
	  //下：初期質量(つまり最小天体の質量)の10倍以上の大きさになったらrunaway bodyとする
    //Information of Jupiter is in a header fout_ps.
	  
	fout_ps << time_year << " " << eJ << " " << aJ << " " << MJ << '\n' << '\n' << '\n' ;
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
	      fout_ps << e[i] << " " << a[i] << " " << m[i]  << '\n';
	//      H[i]=a[i]*pow(m_k[i]/(3.0*1.0),1.0/3.0);
	    }
	    else if(m[i] >= min_mass*10.0 && m[i] != MJ){
//	      H[i]=a[i]*pow(m[i]/(3.0*M_sun),1.0/3.0);
	//      H[i]=a[i]*pow(m_k[i]/(3.0*1.0),1.0/3.0);
	      //             fout_pp << e[i] << " " << a[i]/one_au << " " << m[i] << " " << H[i]/one_au << " " << 5.0*H[i]/one_au << '\n';
	      //4行目はポイントサイズ
	      fout_pp << e[i] << " " << a[i] << " " << m[i] << " " << m[i]/(pow(10.0,24.0)) << " " <<  time_year << '\n';
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
//Property of second mass
//--------------------

	  for(i=0; i<=n-1; i++)
	    {
	    if(m[i] == second)
		{
		  asecond = a[i];
		  esecond = e[i];
		}

	    }

	  
	  mean_mass =( accumulate(&m[0], &m[n-1],0.0) - MJ )/n  ;//mean_mass exclude Jupiter and Sun
//	  cout <<" mean_mass  = " << mean_mass  << '\n'; 
//	  mout << time_year  << "  " << max << "  " << mean_mass  << "  " << second  << "  " << asecond << "  " <<  esecond  << "  " << resonance_angle  << "  " <<  n   << "  " << Npp  <<'\n';//g
	  mout << time_year  << "  " << max << "  " << mean_mass  << "  " << second  << "  " << asecond << "  " <<  esecond   << "  " <<  n   << "  " << Npp  <<'\n';//g



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
   dice(n,value,number,time_year,file_name_cum);
	      
	  for(i=1;i<=n-1;i++)
	    {
		
// Sigma of e and i
	if(m[i] < 1.0e+30){
		sigma_e_rms += pow(e[i],2.0);
		sigma_I_rms += pow(I[i],2.0); 
	
		sigma_e_mean += e[i];
		sigma_I_mean += I[i]; 
		      
//       cout <<"a[7]=" << a[7] << '\n';
//       cout <<"H[7]=" << H[7] << '\n';
//       cout <<"sigma_e_rms_hill =" << sigma_e_rms_hill << '\n';
		sigma_e_rms_hill += pow(e[i]/(H[i]/a[i]),2.0);
		sigma_I_rms_hill += pow(I[i]/(H[i]/a[i]),2.0); 

		n_field = n-1;
			}
      
		
	    }

	  /*
	  //all
	   for(i=0;i<=n-1;i++)
	    {
	     
		  sigma_e += pow(e[i],2.0);
		  sigma_I += pow(I[i],2.0); 
		  sigma_e_mean += e[i];
		  sigma_I_mean += I[i]; 
	    }
	  */
	  e_rms = sqrt(sigma_e_rms/n_field);//root mean square
	  I_rms = sqrt(sigma_I_rms/n_field);
	  e_mean = sigma_e_mean/n_field;//mean
	  I_mean = sigma_I_mean/n_field;
	  e_rms_hill = sqrt(sigma_e_rms_hill/n_field);//root mean square
	  I_rms_hill = sqrt(sigma_I_rms_hill/n_field);

//    cout <<"e_rms_hill=" << sigma_e_rms_hill << '\n';
	  eout << time_year << "  " <<  e_rms << "  " <<  I_rms << "  " <<  e_mean << "  " <<  I_mean << "  " << e_rms_hill << "  " << I_rms_hill << '\n';
	}

	}

    }	
 
      return 0;

    
}


