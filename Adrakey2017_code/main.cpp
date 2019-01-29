#include <iterator>
#include <algorithm>
#include "spatemp.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"
#include "incgammabeta.h"
#include "likelihood.h"
//#include <ctime.h>
#include "Selke.h"
//#include "mcmc_sptemp_new.h"
#include "spatio_temporal_simu.h"
//#include "Simul_remov_ind.h"
//#include "Distri_savewholeregion.h"

using namespace std;


int main(){
Matdoub Coo,M,Selke,time,C,Param,time_interval,Map,selke;
Vecdoub Inf_time,Inf_time_obs,Q,map,threat;
Vecint  Indx_ind,a;
  
 double alpha=.08; double beta=.000007;double epsilon=5e-05; int N=1000; double Tmax=360; double T2=360; int nbur=0000;int thin=1;int nits=10001;int NI=128; double dt=100; int n=(nits-nbur)/thin; double dtype=0.005;
Vecdoub testim{461};Vecdoub l{0.00, 0.075, 0.15,0.225, 0.30,0.375,0.45,0.525,0.6, 0.675, 0.75};
Vecdoub L{0.00, 0.075, 0.15,0.225, 0.30,0.375,0.45,0.525,0.6, 0.675, 0.75}; Matdoub M1;int k;
Vecdoub tim_snapshot{0,30,60,90,120,150,180,210,240,270};
//ofstream ofile("par.txt");
ifstream file;
ofstream ofile;
string ofile1="par_new_pred460.txt";
string ofile2="time_pred460.txt";
string ofile3="selke_pred460.txt";
string simulation="simul54.txt"; 
string file1="Coo.txt";
file_mat(file,N,2,Coo,file1);
/*--------------------------------------     Draw sample from the psoterior ----------*/
n=100001; T2=500;
//   string time1="time_pred1.txt";
//   string par1="par_pred1.txt";
// // // //string selke1="selke_predSI1I2_500.txt";
//   file_mat(file,n-1,N,time,time1);
//    file_mat(file,n-1,5,Param,par1);
// // // // //file_mat(file,n-1,N,selke,selke1);
//   int samp=1000;
//   Matdoub Para(samp,Vecdoub(5));
//  Matdoub Selk(samp,Vecdoub(N));
//   Matdoub Time(samp,Vecdoub(N));
//   Vecdoub inf_time(N);
//  string Para_draw="Param_draw500.txt";
// string Selke_draw="Selked_raw500.txt";
//  string Time_draw="Time_draw500.txt";
// for(int i=0;i<samp;i++){
//   Ran runif9(rand());
// // Ran runif(rand());
// //  int j=1+runif.int64() % (800000-1);
// int j=1+runif9.int64() % (n-1);
// for(int k=0;k<N;k++){
//   inf_time[k]=time[j][k];
// }
// Para[i][0]=Param[j][0];Para[i][1]=Param[j][1];Para[i][2]=Param[j][2];Para[i][3]=Param[j][3];Para[i][4]=Param[j][4];
//  for(int mm=0;mm<N;mm++){
//    Selk[i][mm]=Selke_func(N,Coo,inf_time,Para[i][0],Para[i][1],Para[i][2],T2,mm,dtype);
//   // Time[i][mm]=time[j][mm];
//  }
//  //cout<<time[i][0]<<endl;
// }
// ou_file(ofile,samp,5,Para,Para_draw);
//  ou_file(ofile,samp,N,Selk,Selke_draw);
//ou_file(ofile,samp,N,Time,Time_draw);
/*--------------------------------------     Maps-----------------------------*/
n=100001; T2=500;
// string file7="time_pred1.txt";
// string file8="par_pred1.txt";
// string file9="Selke_draw.txt";
//   string ofile8="Map500.txt";
// //
// file_mat(file,n,N,time,file7);
// file_mat(file,n,5,Param,file8);
//file_mat(file,n,N,selke,file9);
//
//// //cout<<Param[1][0]<<endl;
// Maps(Coo,time,Param,dtype,N,T2,n-1,Map);
// ou_file(ofile,4,N,Map,ofile8);

///*-------------------------------------   Distribution of the number of hosts rescued ------------------------*/
n=1000;
string Para_draw="Param_draw500.txt";
string Selke_draw="Selked_raw500.txt";
// string Time_draw="Time_draw500.txt";
file_mat(file,n,5,Param,Para_draw);
file_mat(file,n,N,selke,Selke_draw);

/*-------------------------------------   Epidemic distribution ------------------------*/
Tmax=500;
int Nepidem2=1000;
Vecint M3(Nepidem2);
Matint M5(Nepidem2,Vecint(4));
Vecdoub testim1{2*Tmax};

#pragma omp parallel for
for(int i=0;i<Nepidem2;i++){
  M3[i]=SI_mimulation_test2(Coo,selke,Param,N,Tmax,dt,{1000},100,dtype,nt,threat,L,l,i);
  //cout<<M3[i]<<endl;
}


string ofile20="M3.txt";
ou_file_vec(ofile,Nepidem2,M3,ofile20);

//file_vec(file,Nepidem2,M3,ofile20);


// varying time
// for(int j=0;j<4;j++){
//   cout<<j<<endl;
// #pragma omp parallel for
// for(int i=0;i<Nepidem2;i++){
//   M5[i][j]=SI_mimulation_test2(Coo,selke,Param,N,Tim[j],dt,testim1,100,dtype,{42,17, 8, 20, 13},threat,L,l,i);
//   // cout<<M3[i]<<endl;
// }
// }
string ofile21="M.txt";
//ou_file(ofile,Nepidem2,4,M5,ofile21);
//file_mat(file,Nepidem2,4,M5,ofile21);


/*--------------------------------------Simulation test -----------------------------*/
string simulation1="simul4.txt";
string episiz="episiz4.txt";
Tmax=360; alpha=.04;  beta=.000007; epsilon=5e-05;
Vecint MM(1000);
string selke1="Q.txt";
file_mat(file,N,1000,Q,selke1);
Matdoub threat1(1000,Vecdoub(4));
// #pragma omp parallel for
//   for(int i=0;i<1000;i++){
//     //cout<<i<<endl;
//     MM[i]=SI_mimulation_test2(Coo,Q,Para,N,Tmax,dt,{1000},nr,dtype,nt,threat,L,l,i);
//     //cout<<MM[i]<<endl;
//   }
//Tmax=500;
// SI_mimulation_test(Coo,Q,Param,N,Tmax,dt,{1000},nr,dtype,nt,threat1,L,l,M1,k,3);
//SI_mimulation_test_one(Coo,Q,alpha,beta,epsilon,N,Tmax,dt,{1000},nr,dtype,nt,M1,k,0);
//int ll=SI_mimulation_test2(Coo,Q,alpha,beta,epsilon,N,Tmax,dt,testim,nr,dtype,nt,threat,L,l);
//int ll=SI_mimulation_test1(Coo,Q,alpha,beta,epsilon,N,Tmax,dt,{58},dtype,0);
//cout<<"ll="<<ll<<endl;
//ou_file(ofile,k,6,M1,simulation1);
//ou_file_vec(ofile,1000,MM,episiz);
/*--------------------------------------  Simulation ----------------------*/
n=1000;
T2=500;
string file10="Map.txt";
file_mat(file,4, N,Map, file10);
string Para_draw="Param_draw.txt";
string Selke_draw="Selke_draw.txt";
string Time_draw="Time_draw.txt";
file_mat(file,n,5,Param,Para_draw);
file_mat(file,n,N,selke,Selke_draw);
int Nepidem=20; int epidem=0;
int NR=50; Matint M2; 
Tmax=500;

threat.resize(N);
Vecdoub selke2(N);
for(int j=0;j<N;j++){
  threat[j]=Map[1][j];
  selke2[j]=selke[epidem][j];
}

T2=500;
Vecint Nind{128,200,300,400,500};
Vecint Tim{510,520,530,540};
vector<vector<string>> v1;
v1.resize(5);
for(int i=0;i<5;i++){
  v1[i].resize(4);
}

// Epidemic size
int Nepidem2=1000;
Vecint M3(Nepidem2);
Matint M5(Nepidem2,Vecint(4));
Vecdoub testim1{2*Tmax};
Vecdoub Selke1(N);
Vecint rank;

string ofile20="M3.txt";;

file_vec(file,Nepidem2,M3,ofile20);



/*------------------Distribution of the hosts saved --------------*/
Nepidem=1000;
T2=500;
Matint Mat(Nepidem,Vecint(5));
Matint Mat1(Nepidem,Vecint(5));
map.resize(N);
for(int k=0;k<N;k++){
  map[k]=Map[0][k];
}


for(int i=0;i<5;i++){
#pragma omp parallel for
  for(int j=0;j<Nepidem;j++){
    Vecint Mi(2);
    SI_mimulation_test4(Coo,selke,Param, N,T2,dt,testim,1,dtype,{Nind[i]},map,{0,0.75},{0,0.75},j,Mi);
    Mat[j][i]=M3[j]-Mi[1];
    Mat1[j][i]=Mi[0];
  }
}
string ofile13="Distup_R.txt";
string ofile14="Distrank_R.txt";
ou_file(ofile,Nepidem,5,Mat1,ofile13);
ou_file(ofile,Nepidem,5,Mat,ofile14);

    return EXIT_SUCCESS;
 
}
