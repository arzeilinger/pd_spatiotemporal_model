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
#include "mcmc_sptemp.h"
#include "spatio_temporal_simu.h"
#include "Simul_remov_ind.h"
#include "Maps.h"
#include "SA.h"
#include "SA_each.h"
// g++ eg1.cpp -std=c++0x
using namespace std;


int main(){
Matdoub Coo,M,Selke,time,C,Param,time_interval,Map,selke;
Vecdoub Inf_time,Inf_time_obs,Q,map,threat;
Vecint  Indx_ind,a;
  
 double alpha=.09; double beta=.0000053;double epsilon=8.2e-05; int N=1000; double Tmax=360; double T2=500; int nbur=0000;int thin=1;int nits=1000001;int NI=128; double dt=0; int n=(nits-nbur)/thin; double dtype=0.005;
Vecdoub testim{361};int nr=6;int nt=7;Vecdoub l{-424.9,-424.5,-424.1,-423.7,-423.3,-422.9};/*{2879.0,2879.2,2879.4,2879.6,2879.8,2880,2880.2,2880.4}*/;Vecdoub L{2879.0,2879.4,2879.8,2880.2,2880.6}/*l{-424.9,-423.7,-424.5,-424.3,-424.1,-423.9,-423.7,-423.5,-423.3,-423.1}*/; Matdoub M1;int k;
Vecdoub tim_snapshot{0,30,60,90,120,150,180,210,240,270};
ifstream file;
//ofstream ofile("par.txt");
ofstream ofile;
string ofile1="par_pred3601.txt";
string ofile2="time_pred3601.txt";
string ofile3="selke_pred3601.txt";
string file1="Coo.txt";
string simulation="simul1.txt"; 
string file2="Inf_tim.txt";
string file3="Indx_ind.txt";
string file4="time_interval.txt";   
string file5="a.txt";
file_mat(file,N,2,Coo,file1);
file_mat(file,N,2,time_interval,file4);
file_vec(file,N,Inf_time, file2);
file_vec(file, N,Indx_ind, file3);
for(int i=0;i<N;i++){
  if(Inf_time[i]>=Tmax){
    Inf_time[i]=T2;
  }
}


double lnew, lold, t_new, pacc;
Vecdoub Inf_tim_new(N);
int Ni=NI,NI_new=0;
int ni=(nits-nbur)/thin; // size after burning (nits-nburn)/nthining o be determined
int acc=1;
Selke.resize(ni);time.resize(ni); Param.resize(ni);
for(int i=0;i<ni;i++){
  Selke[i].resize(N);time[i].resize(N); Param[i].resize(5);  // Vector for Selke, risk and challenge storage
}
Vecdoub Alpha(nits),Beta(nits),Epsilon(nits),Dt(nits); //vector for parameter storage
Beta[0]=7e-06; Alpha[0]=0.1 ;Epsilon[0]=5e-05;Dt[0]=dt; int c=1000; int d=10;//   Initial values
int u=1;
double eps=1.000000;
lold=Likel(Coo,Inf_time,Indx_ind,Alpha[0],Beta[0],Epsilon[0],dtype,N, Tmax, T2,NI);
//--------------------
for(int i=1;i<nits;i++){
  Normaldev rnorm1(Beta[i-1],.00001,rand());
  double Beta_new=rnorm1.dev();
  double pacc=0;
  if(Beta_new>0&& Beta_new<1){
    lnew=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],abs(Beta_new),Epsilon[i-1],dtype,N, Tmax, T2,NI);
    
    pacc=exp(lnew-lold);
    Ran runif(rand());
    //    cout<< NI<<endl;
    if(pacc>runif.doub()){
      Beta[i]=abs(Beta_new);    
      acc+=1;
    }  
    else{
      Beta[i]=Beta[i-1];                                 //update Beta
    }             
  }
  else{
    Beta[i]=Beta[i-1];
  }
  
  //update Epsilon
  
  Normaldev rnorm2(Epsilon[i-1],0.00005,rand());
  double Epsilon_new=rnorm2.dev();
  if(Epsilon_new>0 && Epsilon_new<1){
    lnew=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],Beta[i],abs(Epsilon_new),dtype,N, Tmax, T2,NI);
    // lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],Beta[i],Epsilon[i-1],dtype,N, Tmax, T2,NI);
    pacc=exp(lnew-lold);
    Ran runif7(rand());
    if(pacc>runif7.doub()){
      Epsilon[i]=abs(Epsilon_new);  
      lold=lnew;
    }  
    else{
      Epsilon[i]=Epsilon[i-1];                                
    }             
  }
  else{
    Epsilon[i]=Epsilon[i-1];
  }
  
  //update Alpha   
  
  Normaldev rnorm3(Alpha[i-1],.1,rand());
  double Alpha_new=rnorm3.dev();
  if(Alpha_new>0 && Alpha_new<1){
    lnew=Likel(Coo,Inf_time,Indx_ind,abs(Alpha_new),Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
    //lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
    pacc=exp(lnew-lold);
    Ran runif8(rand());
    if(pacc>runif8.doub()){
      Alpha[i]=abs(Alpha_new);   
      lold=lnew;
    }  
    else{
      Alpha[i]=Alpha[i-1];                                 //update Beta
    }             
  }
  else{
    Alpha[i]=Alpha[i-1];
  }
  Dt[i]=Dt[i-1];
  //Update the infection times
  int k=1;
  while(k<=10){
    Ran runif9(rand());
    int j=1+runif9.int64() % (N-1);  // a randomly choosen individual
    //Symptomatic individual
    int  ll=j-1;
    if((time_interval[ll][1]>2) ){        // Move
      Ran runif11(rand());
      t_new= time_interval[ll][0]+runif11.doub()*(time_interval[ll][1]- time_interval[ll][0]);
      copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
      Inf_tim_new[j-1]=t_new;
      lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
      //lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
      pacc=exp(lnew-lold);
      Ran runif12(rand());
      
      if(pacc>runif12.doub()){
        Inf_time[j-1]=t_new;
        lold=lnew;
      }
    }
    else if(time_interval[ll][0]==1 && time_interval[ll][1]==1){
      Ran runif3(rand());
      if(runif3.doub()<0.5){  // Move infection time
        Ran runif4(rand());
        t_new= Tmax-Dt[i] +runif4.doub()*(T2- Tmax+Dt[i]);// new infection time
        copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
        Inf_tim_new[j-1]=t_new;
        lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
        // lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
        pacc=exp(lnew-lold);
        Ran runif5(rand());
        if(pacc>runif5.doub()){
          Inf_time[j-1]=t_new;
          lold=lnew;
        }                                                                                                                 
      }
      else{  // Delete
        // Inf_tim_new.swap(Inf_time);
        copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
        Inf_tim_new[j-1]=T2;
        NI_new=NI-1;
        lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI_new);
        //lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
        pacc=2.0/(Dt[i]-Tmax+T2)*exp(lnew-lold);
        Ran runif6(rand());
        if(pacc>runif6.doub()){
          Inf_time[j-1]=T2;
          NI=NI_new;
          lold=lnew;
          time_interval[ll][0]=0;
          time_interval[ll][1]=0;
        }
      }
    }
    
    else if(time_interval[ll][0]==0 && time_interval[ll][1]==0) {// Add an infection time
      
      //Inf_tim_new.swap(Inf_time);
      Ran runif1(rand());       
      t_new= Tmax-Dt[i] +runif1.doub()*(T2- Tmax+Dt[i]);// new infection time
      copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
      
      Inf_tim_new[j-1]=t_new;
      NI_new=NI+1;
      lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI_new);
      // lold=Likel(Coo,Inf_time,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
      pacc=.5*(T2-Tmax+Dt[i])*exp(lnew-lold);
      Ran runif2(rand());
      if(pacc>runif2.doub()){
        Inf_time[j-1]=t_new;
        lold=lnew;
        NI=NI_new;
        time_interval[ll][0]=1;
        time_interval[ll][1]=1;
      }
    }
    
    
    k++;
    
  }  
  
  //Selke reconstructiion
  if(i>nbur&& i%thin==0){
    //  cout<< i<<endl;
    for(int j=0;j<N;j++){
      time[u-1][j]=Inf_time[j];
      //   Selke[u-1][j]=Selke_func(N,Coo,Inf_time,Alpha[i], Beta[i],Epsilon[i],T2,j,dtype);
    }
    
    Param[u-1][0]=Alpha[i];
    Param[u-1][1]=Beta[i];
    Param[u-1][2]=Epsilon[i];
    Param[u-1][3]=NI;  
    Param[u-1][4]=Dt[i];    
    u=u+1;
    //cout<<time[u-2][0]<<endl;
  }
  if(i%10000==0){
    //cout<<i<<' '<<u<<endl;
    ou_file(ofile,nits-1,5,Param,ofile1);
    ou_file(ofile,nits-1,N,time,ofile2);
    // ou_file(ofile,nits-1,N,Selke,ofile3);
  }
}


return EXIT_SUCCESS;

}
