#include<iostream>
#include<fstream>
#include<iomanip>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <list>
#include <numeric>



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

int main(int argc, char* argv[])
{

  cout <<"You need to delete dat before starting this analysis." << '\n';
  //単位系　cgs
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


  char filename[20];
  char directoryname[100];
  int i;//ファイル内のループ
  int j;//ファイルのループ
  int k;//ディレクトリのループ

  int K;//使うディレクトリ数

  //  K = 9;
  //  cout << "enter number of directories.\n" ;
  //  cin >> K ;

  K = (int)strtol(argv[1], NULL, 10);

  double time;
  int n; //粒子数
  int n_start; //初期条件での粒子数

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
  double max,second,min,mean_mass;
  double sigma_e = 0.0;
  double sigma_I = 0.0;
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
  f_gas =  strtod(argv[2], NULL);
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
  a_start_au = strtod(argv[3], NULL);
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
        
  
  for(k=1;k<=K;k++)
    {
      sprintf(directoryname,"../../snap%02d",k);
      cout << directoryname << '\n';
	  
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
	      // cout << name << '\n';
	      count++;
	    }   
	}
      closedir(d);
      printf("ファイルが %d 個あります\n", count);
      
      N = count-1;
      //Nは読み込むファイルの数
   
  
      for(j=0;j<=N;j++)
	{

	  sigma_e = 0.0;
	  sigma_I = 0.0; 
	
	  sprintf(filename,"../../snap%02d/snap%05d.dat",k,j);//結果のファイルがどこにあるか確認すること
	  //   cout << filename << '\n';
	  
	  ifstream fin(filename);
	  if(!fin)
	    {
	      cout << "I can't open the file. \n";
	      return 1;
	    }

	  fin >> n >> first_line[0] >>  first_line[1] >> first_line[2] >> first_line[3] >> first_line[4] >> first_line[5] >> first_line[6] >> first_line[7] >> first_line[8] >> first_line[9] >> first_line[10];
	  cout <<"n=" << n << '\n';
	  cout <<"time(year)=" << first_line[0]*(T_unit/one_year)  << '\n';
       
          if( k == 1 && j == 0 )
             {
               n_start = n;
                   }
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
  
	  //その他の情報を読み込んでおくための変数
	  a_x.resize(n);
	  a_y.resize(n);
	  a_z.resize(n);
	  x_1.resize(n);
	  x_2.resize(n);

	  value.resize(n); //出る値
	  number.resize(n);//その値が出た回数

	  n_field = 0;
	  sigma_e = 0.0;
	  sigma_I = 0.0;

  
	  for(i=0;i<=n-1 ;i++)
	    {
	      fin >> id[i] >>m_k[i] >>  x_au[i] >>y_au[i] >> z_au[i] >> v_x_cm[i] >> v_y_cm[i] >> v_z_cm[i] >> a_x[i] >> a_y[i] >> a_z[i] >> x_1[i] >> x_2[i] ;
	      //単位系の換算
	      x[i]=x_au[i]*L_unit;
	      y[i]=y_au[i]*L_unit;
	      z[i]=z_au[i]*L_unit;
    
	      v_x[i]=v_x_cm[i]*V_unit; //*30*pow(10.0,3.0);
	      v_y[i]=v_y_cm[i]*V_unit; //*30*pow(10.0,3.0);
	      v_z[i]=v_z_cm[i]*V_unit; //*30*pow(10.0,3.0);
    
	      m[i]= m_k[i]*M_unit; //*M_sun;
    
	      //軌道要素の計算
	      mu[i]=G_unit*(1.0*M_unit+m[i]);
  
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
	      v_kep = sqrt(G_unit*M_sun/a[i]);  

	   
	    }

	  //inclinationは右向きに傾いているときに負の値をとるので、絶対値をとる
	  for(i=0 ;i<=n-1 ;i++)
	    {
	      I_abs[i] = abs(I[i]);
	    }
	  //relative velocity and max mass
	  max = *max_element(m.begin(), m.end());
	  min = *min_element(m.begin(), m.end());
          // min = *min_element(&m[0], &m[n-1]);

	  //second largest mass
	  for(i=0; i<=n-1; i++)
	    { second = min;
	      if(second < m[i] && m[i] != max)
		{
		  second = m[i];
		}
	    }

	  
	//  mean_mass = accumulate(&m[0], &m[n-1],0.0)/n  ;
	  mean_mass = accumulate(m.begin(), m.end(),0.0)/n;
          cout <<" mean_mass  = " << mean_mass  << '\n'; 
	  mout << first_line[0]*(T_unit/one_year)  << "  " << max << "  " << mean_mass  << "  " << second  << "  " <<  n   <<'\n';//g
	  mout_sigma << first_line[0]*(T_unit/one_year)/pow(a_start_au,3.0/2.0)*Sigma_solid << "  " << max << "  " << mean_mass  << "  " <<  n  << '\n';//g //time[Tkep sigma_solid]

	  cout <<"min=" << min << '\n';
	  // cout <<"n=" << n << '\n';



	  //sort and cumulative mass

	  for(i=0 ;i<=n-1 ;i++)
	    {
	      value[i] = m[i] ;
	      number[i] = 1.0;
	    }
	  stable_sort(value.begin(),value.end(),greater<double>());
	  time = first_line[0]*(T_unit/one_year) ;
	  dice(n,value,number,time,a_start_au);

	      
	  for(i=0;i<=n-1;i++)
	    {
	      if( k == 1 && j == 0 || n == n_start)
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

		}
	      else
		{
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
		      sigma_e += pow(e[i],2.0);
		      sigma_I += pow(I[i],2.0); 
	
		      n_field += 1;
		      
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
	  Sigma_gas = 2400.0 * f_gas *pow(a_start_au,-3.0/2.0)*exp(-1.0*first_line[0]*(T_unit/one_year)/tau_gas) ;
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

    
	  eout << first_line[0]*(T_unit/one_year) << "  " <<  e_m[k] <<  '\n';
	  fout << first_line[0]*(T_unit/one_year) << "  " << f << '\n';
	  e_analytical_out << first_line[0]*(T_unit/one_year)<< "  " << de2_cd_1 << "  "<< de2_cd_2 << "  " << de2_gas << "  " << de2_vs_1 << "  " << de2_vs_2 << "  " << de2_df <<  '\n';
	}



    }	
 
      return 0;

    
}

