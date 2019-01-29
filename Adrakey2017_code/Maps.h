//function to create maps i.e measure use for the control

void Maps(Matdoub Coo, Matdoub time,Matdoub paramet, int dtype, int N, double T2,int ni,Matdoub &map){
  Matdoub Risk(ni,Vecdoub(N));
  Matdoub chal(ni,Vecdoub(N));
  double s,d,s1,s2,s3;
//#pragma omp parallel for
  
  for(int i=0;i<ni;i++){
   // cout<<i<<endl;
    for(int j=0;j<N;j++){
    
      if(time[i][j]<T2){ // Host j infected
        Risk[i][j]=1;
      }
      else{
        Risk[i][j]=0;
        
      }
 //     cout<<Risk[i][j]<<endl;
      s=0;
      for(int k=0;k<N;k++){
          if(k!=j && time[i][k]>=T2){
            d= sqrt(pow(Coo[j][0]-Coo[k][0],2)+ pow(Coo[j][1]-Coo[k][1],2));
            s+=paramet[i][1]*kernel(d,paramet[i][0],dtype);
          }

        }
      chal[i][j]=s;
      
    }
  }
  
map.resize(4);
for(int i=0;i<4;i++){
  map[i].resize(N);
}

for(int i=0;i<N;i++){
  s=0;s1=0;s2=0;
  for(int j=0;j<ni;j++){
    s+=Risk[j][i];
    s1+=chal[j][i];
    s2+=Risk[j][i]*chal[j][i];
  }
  map[0][i]=s/ni;                     // Risk
  map[1][i]=s1/ni;                    //  Hazard
  map[2][i]=s2/ni;                   // Threat
  map[3][i]=map[0][i]*map[1][i];
  
}
}
