double Selke_func(int N, Matdoub Coo,Vecdoub Inf_time,double alpha, double beta,double epsilon, double T2,int j,double dtype){
double selke;
double d;
selke= epsilon*Inf_time[j];
for(int k=0; k<N;k++){
  if(Inf_time[j]>Inf_time[k]){
    d= sqrt(pow(Coo[j][0]-Coo[k][0],2)+ pow(Coo[j][1]-Coo[k][1],2));
    selke= selke + beta*kernel(d,alpha,dtype)*(Inf_time[j]-Inf_time[k]);
  }
}
 
if(Inf_time[j]>=T2) {
 Ran runif(rand());
 selke= selke + runif.doub();
}
 
 return selke;
}
