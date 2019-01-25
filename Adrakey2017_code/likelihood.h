inline double Likel(Matdoub Coo, Vecdoub Inf_time, Vecint Indx_ind,double alpha, double beta, double epsilon, double dtype, int N, double Tmax, double T2, int NI ){
 double Suc_pres,Inf_pres,Inf_pres1=0; // initialization
 // Here we assume that the loglikelihoo is in the form Inf_pres-Inf_pres -Suc_pres
 double S,P,m=0; 
 Matdoub dat(NI,Vecdoub (4));
// Inf_pres=epsilon*dat[0][0];    // First infection
 data(Coo, Inf_time, Indx_ind,dat,NI); // data rearangement by colum: time, x,y and index resp

 double d;
 //cout<<S<<endl;
 for(int k=0;k<=N;k++){
// cout<<S<<endl;
  if(k==N){
    for(int j=0;j<N;j++){
       S+= epsilon*Inf_time[j];
  }
  }
  else if(Inf_time[k]<T2){
//#pragma omp parallel for reduction(+:S)
    for(int j=0;j<N;j++){
       if(Inf_time[j]>Inf_time[k]){
         d= sqrt(pow(Coo[j][0]-Coo[k][0],2)+ pow(Coo[j][1]-Coo[k][1],2));
         S=S+beta*kernel(d,alpha,dtype)*(Inf_time[j]-Inf_time[k]);
       }
   }
   
   }
  }
//   cout<<S<<endl;
 for(int k=1;k<NI;k++){
     m=0;
     for(int j=0;j<k;j++){
      d= sqrt(pow(dat[j][1]-dat[k][1],2)+ pow(dat[j][2]-dat[k][2],2));      
      m+=beta*kernel(d,alpha,dtype) ;  
 //     cout<<d<<endl;
     }
     P+=log(m+epsilon);
 }
//cout<<P<<endl;
 return(P-S);
 

}
