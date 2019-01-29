void SII(int nbur,int thin,int nits,Matdoub Coo,Vecdoub Inf_time, Vecint Indx_ind,double dtype, int N, double Tmax, double T2,int NI,Matdoub &Selke,Matdoub &time,Matdoub &Param,double dt,Matdoub time_interval)
{

// Define some ..
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
Beta[0]=6.1479804059e-06; Alpha[0]=0.0771343098 ;Epsilon[0]=9.26915e-05;Dt[0]=dt; int c=1000; int d=10;//   Initial values
int u=1;
double eps=1.000000;
lold=Likel(Coo,Inf_time,Indx_ind,Alpha[0],Beta[0],Epsilon[0],dtype,N, Tmax, T2,NI);
//--------------------
for(int i=1;i<nits;i++){
  Normaldev rnorm1(Beta[i-1],.00001,rand());
  double Beta_new=rnorm1.dev();
  double pacc=0;
  if(Beta_new>0 && Beta_new<1){
     lnew=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],Beta_new,Epsilon[i-1],dtype,N, Tmax, T2,NI);

    pacc=exp(lnew-lold);
    Ran runif(rand());
//    cout<< NI<<endl;
    if(pacc>runif.doub()){
     Beta[i]=Beta_new;
      lold=lnew;
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
    lnew=Likel(Coo,Inf_time,Indx_ind,Alpha[i-1],Beta[i],Epsilon_new,dtype,N, Tmax, T2,NI);
    double pacc=exp(lnew-lold);
    Ran runif7(rand());
    if(pacc>runif7.doub()){
     Epsilon[i]=Epsilon_new;
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
   lnew=Likel(Coo,Inf_time,Indx_ind,Alpha_new,Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
    double pacc=exp(lnew-lold);
    Ran runif8(rand());
    if(pacc>runif8.doub()){
     Alpha[i]=Alpha_new;
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
//}
  //Update the infection times
  int k=1;
  while(k<1){
    Ran runif9(rand());
   int j1=1+runif9.int64() % (NI-1);  // a randomly choosen individual
   int j=Indx_ind[j1];
    //Symptomatic individual
 int  ll=j-1;
   if((time_interval[ll][1]>1) ){        // Move
  //  cout<< ll<<endl;
        Ran runif11(rand());
       t_new= time_interval[ll][0]+runif11.doub()*(time_interval[ll][1]- time_interval[ll][0]);
      // cout<<time_interval[ll][1]<<" "<<j<<" ";
       copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
       //Inf_tim_new.swap(Inf_time);
       Inf_tim_new[j-1]=t_new;
       lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
       pacc=exp(lnew-lold);
       Ran runif12(rand());

      if(pacc>runif12.doub()){
         Inf_time[j-1]=t_new;
       }
      }
   else if(time_interval[ll][0]==1 && time_interval[ll][1]==1){
 //  cout<< ll<<" "<<k<<endl;
         Ran runif3(rand());
      if(runif3.doub()<0.5){  // Move infection time
         Ran runif4(rand());
         t_new= Tmax-Dt[i] +runif4.doub()*(T2- Tmax+Dt[i]);// new infection time
         copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
        // Inf_tim_new.swap(Inf_time);
         Inf_tim_new[j-1]=t_new;
         lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI);
         pacc=exp(lnew-lold);
          Ran runif5(rand());
        if(pacc>runif5.doub()){
          Inf_time[j-1]=t_new;
         }
      }
      else{  // Delete
       // Inf_tim_new.swap(Inf_time);
       copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
        Inf_tim_new[j-1]=T2;
        NI_new=NI-1;
         lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI_new);
         pacc=2.0/(Dt[i]-Tmax+T2)*exp(lnew-lold);
         Ran runif6(rand());
         if(pacc>runif6.doub()){
           Inf_time[j-1]=T2;
           NI=NI_new;
           time_interval[ll][0]=0;
           time_interval[ll][1]=0;
         }
      }
     }

    else if(time_interval[ll][0]==0 && time_interval[ll][1]==0) {// Add an infection time

       //Inf_tim_new.swap(Inf_time);
       Ran runif1(rand());
       t_new= Tmax-Dt[i] +runif1.doub()*(T2- Tmax+Dt[i]);// new infection time
      // Inf_tim_new.swap(Inf_time);
      copy(Inf_time.begin(), Inf_time.end(), Inf_tim_new.begin());
   //    cout<< ll<<" ";

       Inf_tim_new[j-1]=t_new;
       NI_new=NI+1;
       lnew=Likel(Coo,Inf_tim_new,Indx_ind,Alpha[i],Beta[i],Epsilon[i],dtype,N, Tmax, T2,NI_new);
       pacc=.5*(T2-Tmax+Dt[i])*exp(lnew-lold);
    //   cout<<(T2-Tmax+dt)*exp(lnew-lold)<<" "<<pacc<<endl;
       Ran runif2(rand());
       if(pacc>runif2.doub()){
         Inf_time[j-1]=t_new;
         NI=NI_new;
         time_interval[ll][0]=1;
         time_interval[ll][1]=1;
         lold=lnew;
       }
     }


    k++;

  }

  //Selke reconstructiion
  //cout<<s<<endl;
  //  cout<< i<<endl;
    //for(int j=0;j<N;j++){
     // time[i][j]=Inf_time[j];
      //Selke[u-1][j]=Selke_func(N,Coo,Inf_time,Alpha[i], Beta[i],Epsilon[i],T2,j,dtype);
    //}

  Param[i][0]=Alpha[i];
  Param[i][1]=Beta[i];
  Param[i][2]=Epsilon[i];
  Param[i][3]=NI;
  Param[i][4]=Dt[i];
 // cout<<i<<" "<<NI<<endl;
  
}

//cout<<acc<<endl;
}
