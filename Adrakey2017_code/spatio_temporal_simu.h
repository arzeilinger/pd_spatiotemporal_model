//---------------- This functions are intended to simulate an epidemic SI process with cryptic hosts amond the population


//Identiy individual in all region

void Ind_region(Vecint X,Matdoub Coo, Vecdoub L,Vecdoub l,Matint &Ind_test, int N){
  int n=L.size();
  int m=l.size();
  int k=0;int u;
  
  Ind_test.resize(N);
  for(int i=0;i<N;i++){
    Ind_test[i].resize((n-1)*(m-1));
  }
  
  for(int i=0;i<(n-1);i++){
    for(int j=0;j<(m-1);j++){
      u=0;
      for(int r=0;r<X.size();r++){
        //      if(i==1 && j==1){
        //       cout<<Coo[X[r]-1][0]<<" "<<l[j]<<" "<<Coo[X[r]-1][1]<<" "<<l[j+1]<<" "<<L[i]<<" "<<L[i+1]<<endl;
        //      }
        
        if((Coo[X[r]-1][0]>=l[j] && Coo[X[r]-1][0]<l[j+1])&&(Coo[X[r]-1][1]>=L[i] && Coo[X[r]-1][1]<L[i+1])){
          Ind_test[u][k]=X[r];
          u+=1;
        }
        
      }
      k+=1;
    }
  }
}


// Fucntion to simulate an epidemic with removal due to control

void SI_mimulation_test(Matdoub Coo,Matdoub Selke,Matdoub para, int N, double Tmax,  double dt, Vecdoub testim, int nr,int dtype,int nt,Matdoub Map,Vecdoub L,Vecdoub l,Matdoub &M,int &k,int m1){
  //Initial data
  Vecdoub selke(N); Vecdoub map(N);
  double alpha=para[m1][0]; double beta=para[m1][1]; double epsilon=para[m1][2];
  for(int i=0;i<N;i++){
    selke[i]=Selke[m1][i];
    map[i]=Map[2][i];
  }
  double t=0; double s, d,Min1,tt;
  Vecdoub::iterator Min;
  Matdoub M1(2*N,Vecdoub(6));Matdoub M2(2*N,Vecdoub(5));
  Matint  inf_indx(N,Vecint(2));//time of infection and indx of individual infected
  Vecint suc_indx(N); Vecint Rem_indx(N);Vecint X(N);
  Vecdoub Dt(N); Vecdoub Risk(N);  //Pressure
  int n=L.size();
  int m=l.size();
  Vecdoub threat((n-1)*(m-1)); Vecdoub threat_cop((n-1)*(m-1)); Vecdoub reg_map(N);
  for(int i=0;i<N;i++){
    suc_indx[i]=X[i]=i+1;
  }
  int r=0;int rm=0;
  int indx,v,ll;
  Matint Ind_test;
  int u;
  k=0;
  
  while((t<=Tmax)){
    // cout<<t<<" "<<k<<" "<<testim[0]<<endl;
    if(k==0){ // first infection
      for(int i=0;i<N;i++){
        Dt[i]=selke[i]/epsilon;
        Risk[i]=epsilon;
        
      }
      Min=min_element(Dt.begin(),Dt.end());
      indx=distance(Dt.begin(),Min);
      //   u=suc_indx[0];
      //  cout<<Dt[indx]<<endl;
      for(int i=0;i<testim.size();i++){
        testim[i]=testim[i]-Dt[indx];
      }
      t+=Dt[indx];
      Min1=Dt[indx];
      v=0;
    }
    else{
      
      if(r==0){
        for(int i=0; i<N;i++){
          if(suc_indx[i]>=1){ 
            Risk[i]=epsilon;
            Dt[i]=selke[i]/epsilon;
          }
          else{
            Dt[i]=Tmax;
          }
          
        }
        v=0;
      }
      else{
        for(int i=0;i<N;i++){
          if(k==12){
          }
          
          
          if(suc_indx[i]>=1){ // susceptibles
            
            s=0;
            for(int j=0;j<N;j++){
              if(inf_indx[j][1]!=0){ // infected 
                d=sqrt(pow(Coo[i][0]-Coo[j][0],2)+ pow(Coo[i][1]-Coo[j][1],2));
                s=s+beta*kernel(d,alpha,dtype);
              }
              
            }
            Risk[i]=s+epsilon;
            Dt[i]=selke[i]/(s+epsilon);
          }
          else{
            Dt[i]=Tmax;
          }
        }
      }
      Min=min_element(Dt.begin(),Dt.end()); // Min time of infection
      indx=distance(Dt.begin(),Min);
      
      if(testim[0]==0){            
        Min1=Dt[indx];
        t+=Dt[indx];
        v=0;
      }
      else{ // Test or infection
        //   cout<<k<<endl;
        Min1=min(Dt[indx],testim[0]);
        if(Min1==Dt[indx]){ //infection
          t+=Dt[indx];
          for(int i=0;i<testim.size();i++){
            testim[i]=testim[i]-Dt[indx];
          }
          
          v=0;
        }
        else{  //test
          //      cout<<selke[148]<<endl;
          t+=Min1;
          if(testim.size()==1){
            testim[0]=0;
          }
          else{
            testim.erase(testim.begin());
          }
          
          Ind_region(X,Coo,L,l,Ind_test,N);
          
          //          for(int  i=0;i<(n-1)*(m-1);i++){
          //          cout<<i<<"=>";
          //            for(int j=0;j<N;j++){
          //              cout<<Ind_test[j][i]<<" ";
          //            }
          //            cout<<endl;
          //          }
          for(int i=0;i<(n-1)*(m-1);i++){
            s=0;
            tt=0;
            for(int j=0;j<N;j++){
              if(Ind_test[j][i]==0){
                break;
              }
              else{
                s+=map[Ind_test[j][i]-1];
                tt+=1;
              }
            }
            if(tt==0){
              threat[i]=0;
            }
            else{
              threat[i]=s/tt;
            }
            
            //            cout<<threat[i]<<" "<<i<<endl;
            
          }
          
          copy(threat.begin(), threat.end(), threat_cop.begin());
          sort(threat_cop.begin(),threat_cop.end() );            //threat of each region
          //          for(int i=0;i<(n-1)*(m-1);i++){
          //            cout<<threat[i]<<endl;
          //          }
          for(int i=0;i<nr;i++){
            //          cout<<i<<endl;
            if(threat_cop[(n-1)*(m-1)-i-1]>0){                             //presence of host 
              ll=find(threat.begin(),threat.end(),threat_cop[(n-1)*(m-1)-i-1])-threat.begin();
              Vecdoub threat_reves(N); Vecdoub threat_resrt(N);                           // Threat of region to visit
              //  cout<<ll<<endl;
              threat[ll]=0;
              for(int j=0;j<N;j++){
                if(Ind_test[j][ll]==0){   
                  
                  break;
                }
                else{
                  threat_reves[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                  threat_resrt[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                }
                //   cout<<ll<<endl;
              }
              //              int k2=0;
              //              for(int k1=0;k1<N;k1++){
              //              if(threat_resrt[k1]>0){
              //                k2++;
              //              }
              //              }
              //              cout<<k2<<endl;
              sort(threat_resrt.begin(),threat_resrt.end());
              for(int l=0;l<nt;l++){
                //           cout<<threat_resrt[N-l-1]<<endl;
                int ii=find(threat_reves.begin(),threat_reves.end(),threat_resrt[N-l-1])-threat_reves.begin();
                if(threat_resrt[N-l-1]==0){  //susceptibles
                  cout<<N<<endl;
                  break;
                }
                else if(inf_indx[ii][0]!=0&&inf_indx[ii][1]!=0){
                  //Infectives
                  
                  // if((inf_indx[ii][0]<=t)){    //symptomatic
                  inf_indx[ii][0]=0;
                  inf_indx[ii][1]=0;
                  threat_reves[ii]=0;
                  rm+=1;
                  r=r-1;
                  //  }
                }
              }
              
            }
            else{
              break;
            }
            
            
          }
          
          M1[k][0]=t;
          M1[k][1]=0;
          M1[k][2]=0;
          M1[k][3]=0;
          M1[k][4]=r;
          M1[k][5]=rm;
          v=1;
          k+=1;
          //          cout<<k<<" "<<selke[148]<<endl;
        }
      }
      //      }
    }
    if(v==0){
      //        cout<<k<<" "<<indx<<endl;
      
      if(t<=Tmax){
        r+=1;
        M1[k][0]=t;
        M1[k][1]=Coo[indx][0];
        M1[k][2]=Coo[indx][1];
        M1[k][3]=indx+1;
        M1[k][4]=r;
        M1[k][5]=rm;
        k+=1;    
      }      
      inf_indx[indx][1]=indx+1;
      inf_indx[indx][0]=t;
      suc_indx[indx]=0;
      selke[indx]=0;
      //   cout<<r<<endl;
      //Update the thresholds
      
    }
    
    for(int i=0;i<N;i++){
      if(suc_indx[i]!=0){ // susceptibles
        selke[i]=selke[i]-Min1*Risk[i];
      }
    }
    
    
  }
  M.resize(k);
  for(int i=0;i<(k);i++){
    M[i].resize(6);
  }
  
  for(int i=0;i<k;i++){
    
    M[i][0]=M1[i][0];
    M[i][1]=M1[i][1];
    M[i][2]=M1[i][2];
    M[i][3]=M1[i][3];
    M[i][4]=M1[i][4];
    M[i][5]=M1[i][5];
  }
}

// Function to find the size of the population rescued

int SI_mimulation_test2(Matdoub Coo,Matdoub Selke,Matdoub Par, int N, double Tmax,  double dt, Vecdoub testim, int nr,int dtype,int nt,Vecdoub map,Vecdoub L,Vecdoub l, int kk){
  //Initial data
  double alpha=.02;  double beta=.000007; double epsilon=5e-05;
  // double alpha=Par[kk][0];
  // double beta=Par[kk][1];
  // double epsilon=Par[kk][2];
  Vecdoub selke(N);
  for(int i=0;i<N;i++){
    selke[i]=Selke[kk][i];
    // cout<<selke[i]<<' ';
  }
 // cout<<endl;
  
  double t=0; double s, d,Min1,tt;
  Vecdoub::iterator Min;
  Matdoub M1(2*N,Vecdoub(6));Matdoub M2(2*N,Vecdoub(5));
  Matint  inf_indx(N,Vecint(2));//time of infection and indx of individual infected
  Vecint suc_indx(N); Vecint Rem_indx(N);Vecint X(N);
  Vecdoub Dt(N); Vecdoub Risk(N);  //Pressure
  int n=L.size();
  int m=l.size();
  Vecdoub threat((n-1)*(m-1)); Vecdoub threat_cop((n-1)*(m-1)); Vecdoub reg_map(N);
  for(int i=0;i<N;i++){
    suc_indx[i]=X[i]=i+1;
  }
  int r=0;int rm=0; int k=0;
  int indx,v,ll;
  Matint Ind_test;
  int u;
  k=0;
  
  while((t<=Tmax)){
    // cout<<t<<" "<<k<<" "<<testim[0]<<endl;
    if(k==0){ // first infection
      for(int i=0;i<N;i++){
        Dt[i]=selke[i]/epsilon;
        Risk[i]=epsilon;
        
      }
      Min=min_element(Dt.begin(),Dt.end());
      indx=distance(Dt.begin(),Min);
      //   u=suc_indx[0];
      //  cout<<Dt[indx]<<endl;
      for(int i=0;i<testim.size();i++){
        testim[i]=testim[i]-Dt[indx];
      }
      t+=Dt[indx];
      Min1=Dt[indx];
      v=0;
    }
    else{
      
      if(r==0){
        for(int i=0; i<N;i++){
          if(suc_indx[i]>=1){ 
            Risk[i]=epsilon;
            Dt[i]=selke[i]/epsilon;
          }
          else{
            Dt[i]=Tmax;
          }
          
        }
        v=0;
      }
      else{
        for(int i=0;i<N;i++){
          if(k==12){
          }
          
          
          if(suc_indx[i]>=1){ // susceptibles
            
            s=0;
            for(int j=0;j<N;j++){
              if(inf_indx[j][1]!=0){ // infected 
                d=sqrt(pow(Coo[i][0]-Coo[j][0],2)+ pow(Coo[i][1]-Coo[j][1],2));
                s=s+beta*kernel(d,alpha,dtype);
              }
              
            }
            Risk[i]=s+epsilon;
            Dt[i]=selke[i]/(s+epsilon);
          }
          else{
            Dt[i]=Tmax;
          }
        }
      }
      Min=min_element(Dt.begin(),Dt.end()); // Min time of infection
      indx=distance(Dt.begin(),Min);
      
      if(testim[0]==0){            
        Min1=Dt[indx];
        t+=Dt[indx];
        v=0;
      }
      else{ // Test or infection
        //   cout<<k<<endl;
        Min1=min(Dt[indx],testim[0]);
        if(Min1==Dt[indx]){ //infection
          t+=Dt[indx];
          for(int i=0;i<testim.size();i++){
            testim[i]=testim[i]-Dt[indx];
          }
          
          v=0;
        }
        else{  //test
          //      cout<<selke[148]<<endl;
          t+=Min1;
          if(testim.size()==1){
            testim[0]=0;
          }
          else{
            testim.erase(testim.begin());
          }
          
          Ind_region(X,Coo,L,l,Ind_test,N);
          
          //          for(int  i=0;i<(n-1)*(m-1);i++){
          //          cout<<i<<"=>";
          //            for(int j=0;j<N;j++){
          //              cout<<Ind_test[j][i]<<" ";
          //            }
          //            cout<<endl;
          //          }
          for(int i=0;i<(n-1)*(m-1);i++){
            s=0;
            tt=0;
            for(int j=0;j<N;j++){
              if(Ind_test[j][i]==0){
                break;
              }
              else{
                s+=map[Ind_test[j][i]-1];
                tt+=1;
              }
            }
            if(tt==0){
              threat[i]=0;
            }
            else{
              threat[i]=s/tt;
            }
            
            //            cout<<threat[i]<<" "<<i<<endl;
            
          }
          
          copy(threat.begin(), threat.end(), threat_cop.begin());
          sort(threat_cop.begin(),threat_cop.end() );            //threat of each region
          //          for(int i=0;i<(n-1)*(m-1);i++){
          //            cout<<threat[i]<<endl;
          //          }
          for(int i=0;i<nr;i++){
            //          cout<<i<<endl;
            if(threat_cop[(n-1)*(m-1)-i-1]>0){                             //presence of host 
              ll=find(threat.begin(),threat.end(),threat_cop[(n-1)*(m-1)-i-1])-threat.begin();
              Vecdoub threat_reves(N); Vecdoub threat_resrt(N);                           // Threat of region to visit
              //  cout<<ll<<endl;
              threat[ll]=0;
              for(int j=0;j<N;j++){
                if(Ind_test[j][ll]==0){   
                  
                  break;
                }
                else{
                  threat_reves[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                  threat_resrt[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                }
                //   cout<<ll<<endl;
              }
              //              int k2=0;
              //              for(int k1=0;k1<N;k1++){
              //              if(threat_resrt[k1]>0){
              //                k2++;
              //              }
              //              }
              //              cout<<k2<<endl;
              sort(threat_resrt.begin(),threat_resrt.end());
              for(int l=0;l<nt;l++){
                //           cout<<threat_resrt[N-l-1]<<endl;
                int ii=find(threat_reves.begin(),threat_reves.end(),threat_resrt[N-l-1])-threat_reves.begin();
                if(threat_resrt[N-l-1]==0){  //susceptibles
                  //            cout<<N<<endl;
                  break;
                }
                else if(inf_indx[ii][0]!=0&&inf_indx[ii][1]!=0){
                  //Infectives
                  
                  // if((inf_indx[ii][0]<=t)){    //symptomatic
                  inf_indx[ii][0]=0;
                  inf_indx[ii][1]=0;
                  threat_reves[ii]=0;
                  rm+=1;
                  r=r-1;
                  //  }
                }
              }
              
            }
            else{
              break;
            }
            
            
          }
          
          M1[k][0]=t;
          M1[k][1]=0;
          M1[k][2]=0;
          M1[k][3]=0;
          M1[k][4]=r;
          M1[k][5]=rm;
          v=1;
          k+=1;
          //          cout<<k<<" "<<selke[148]<<endl;
        }
      }
      //      }
    }
    if(v==0){
      //        cout<<k<<" "<<indx<<endl;
      
      if(t<=Tmax){
        r+=1;
        M1[k][0]=t;
        M1[k][1]=Coo[indx][0];
        M1[k][2]=Coo[indx][1];
        M1[k][3]=indx+1;
        M1[k][4]=r;
        M1[k][5]=rm;
        k+=1;    
      }      
      inf_indx[indx][1]=indx+1;
      inf_indx[indx][0]=t;
      suc_indx[indx]=0;
      selke[indx]=0;
      //   cout<<r<<endl;
      //Update the thresholds
      
    }
    
    for(int i=0;i<N;i++){
      if(suc_indx[i]!=0){ // susceptibles
        selke[i]=selke[i]-Min1*Risk[i];
      }
    }
    
    
  }
  // cout<<k<<endl;
  //cout<<alpha<<" "<<beta<<" "<<epsilon<<endl;
  return (M1[k-1][4]+M1[k-1][5]);
}

// Function to simulate on epidemic process

void SI_mimulation_testone(Matdoub Coo,Vecdoub selke,double alpha, double beta,double epsilon, int N, double Tmax,  double dt, Vecdoub testim, int nr,int dtype,int nt,Vecdoub map,Vecdoub L,Vecdoub l,Matdoub &M,int &k){
  //Initial data
  double t=0; double s, d,Min1,tt;
  Vecdoub::iterator Min;
  Matdoub M1(2*N,Vecdoub(6));Matdoub M2(2*N,Vecdoub(5));
  Matint  inf_indx(N,Vecint(2));//time of infection and indx of individual infected
  Vecint suc_indx(N); Vecint Rem_indx(N);Vecint X(N);
  Vecdoub Dt(N); Vecdoub Risk(N);  //Pressure
  int n=L.size();
  int m=l.size();
  Vecdoub threat((n-1)*(m-1)); Vecdoub threat_cop((n-1)*(m-1)); Vecdoub reg_map(N);
  for(int i=0;i<N;i++){
    suc_indx[i]=X[i]=i+1;
  }
  int r=0;int rm=0;
  int indx,v,ll;
  Matint Ind_test;
  int u;
  k=0;
  
  while((t<=Tmax)){
    //cout<<t<<" "<<k<<endl;
    if(k==0){ // first infection
      for(int i=0;i<N;i++){
        Dt[i]=selke[i]/epsilon;
        Risk[i]=epsilon;
        
      }
      Min=min_element(Dt.begin(),Dt.end());
      indx=distance(Dt.begin(),Min);
      //   u=suc_indx[0];
      //  cout<<Dt[indx]<<endl;
      for(int i=0;i<testim.size();i++){
        testim[i]=testim[i]-Dt[indx];
      }
      t+=Dt[indx];
      Min1=Dt[indx];
      v=0;
    }
    else{
      if(r==0){
        for(int i=0; i<N;i++){
          if(suc_indx[i]>=1){ 
            Risk[i]=epsilon;
            Dt[i]=selke[i]/epsilon;
          }
          else{
            Dt[i]=Tmax;
          }
          
        }
        v=0;
      }
      else{
        for(int i=0;i<N;i++){
          if(k==12){
          }
          
          
          if(suc_indx[i]>=1){ // susceptibles
            
            s=0;
            for(int j=0;j<N;j++){
              if(inf_indx[j][1]!=0){ // infected 
                d=sqrt(pow(Coo[i][0]-Coo[j][0],2)+ pow(Coo[i][1]-Coo[j][1],2));
                s=s+beta*kernel(d,alpha,dtype);
              }
              
            }
            Risk[i]=s+epsilon;
            Dt[i]=selke[i]/(s+epsilon);
          }
          else{
            Dt[i]=Tmax;
          }
        }
      }
      Min=min_element(Dt.begin(),Dt.end()); // Min time of infection
      indx=distance(Dt.begin(),Min);
      if(k==1){
      }
      if(testim.size()==0){            
        Min1=Dt[indx];
        t+=Dt[indx];
        // cout<<t<<endl;
        v=0;
      }
      else{ // Test or infection
        Min1=min(Dt[indx],testim[0]);
        if(Min1==Dt[indx]){ //infection
          t+=Dt[indx];
          for(int i=0;i<testim.size();i++){
            testim[i]=testim[i]-Dt[indx];
          }
          
          v=0;
        }
        else{  //test
          t+=Min1;
          testim.erase(testim.begin());
          Ind_region(X,Coo,L,l,Ind_test,N);
          
          //          for(int  i=0;i<(n-1)*(m-1);i++){
          //          cout<<i<<"=>";
          //            for(int j=0;j<N;j++){
          //              cout<<Ind_test[j][i]<<" ";
          //            }
          //            cout<<endl;
          //          }
          for(int i=0;i<(n-1)*(m-1);i++){
            s=0;
            tt=0;
            for(int j=0;j<N;j++){
              if(Ind_test[j][i]==0){
                break;
              }
              else{
                s+=map[Ind_test[j][i]-1];
                tt+=1;
              }
            }
            if(tt==0){
              threat[i]=0;
            }
            else{
              threat[i]=s/tt;
            }
            
            //            cout<<threat[i]<<" "<<i<<endl;
            
          }
          
          copy(threat.begin(), threat.end(), threat_cop.begin());
          sort(threat_cop.begin(),threat_cop.end() );            //threat of each region
          //          for(int i=0;i<(n-1)*(m-1);i++){
          //            cout<<threat[i]<<endl;
          //          }
          for(int i=0;i<nr;i++){
            //       cout<<i<<endl;
            if(threat_cop[(n-1)*(m-1)-i-1]>0){                             //presence of host 
              ll=find(threat.begin(),threat.end(),threat_cop[(n-1)*(m-1)-i-1])-threat.begin();
              Vecdoub threat_reves(N); Vecdoub threat_resrt(N);                           // Threat of region to visit
              //     cout<<ll<<endl;
              threat[ll]=0;
              for(int j=0;j<N;j++){
                if(Ind_test[j][ll]==0){   
                  
                  break;
                }
                else{
                  threat_reves[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                  threat_resrt[Ind_test[j][ll]-1]=map[Ind_test[j][ll]-1];
                }
                //   cout<<ll<<endl;
              }
              
              sort(threat_resrt.begin(),threat_resrt.end());
              for(int l=0;l<nt;l++){
                //              cout<<l<<endl;
                if(threat_resrt[N-l-1]==0){  //susceptibles
                  break;
                }
                else{     //Infectives
                  int ii=find(threat_reves.begin(),threat_reves.end(),threat_resrt[N-l-1])-threat_reves.begin();
                  // if((inf_indx[ii][0]<=t)){    //symptomatic
                  inf_indx[ii][0]=0;
                  inf_indx[ii][1]=0;
                  threat_reves[ii]=0;
                  rm+=1;
                  r=r-1;
                  //  }
                }
              }
              
            }
            else{
              break;
            }
            
            
          }
          
          M1[k][0]=t;
          M1[k][1]=0;
          M1[k][2]=0;
          M1[k][3]=0;
          M1[k][4]=r;
          M1[k][5]=rm;
          v=1;
          k+=1;
          
        }
      }
    }
    if(v==0){
      //       cout<<k<<" "<<indx<<endl;
      
      if(t<=Tmax){
        r+=1;
        M1[k][0]=t;
        M1[k][1]=Coo[indx][0];
        M1[k][2]=Coo[indx][1];
        M1[k][3]=indx+1;
        M1[k][4]=r;
        M1[k][5]=rm;
        k+=1;    
      }      
      inf_indx[indx][1]=indx+1;
      inf_indx[indx][0]=t;
      suc_indx[indx]=0;
      selke[indx]=0;
      //   cout<<r<<endl;
      //Update the thresholds
      
    }
    for(int i=0;i<N;i++){
      if(suc_indx[i]!=0){ // susceptibles
        selke[i]=selke[i]-Min1*Risk[i];
      }
    }
    
  }
  M.resize(k);
  for(int i=0;i<(k);i++){
    M[i].resize(6);
  }
  
  for(int i=0;i<k;i++){
    
    M[i][0]=M1[i][0];
    M[i][1]=M1[i][1];
    M[i][2]=M1[i][2];
    M[i][3]=M1[i][3];
    M[i][4]=M1[i][4];
    M[i][5]=M1[i][5];
  }
}


/*----------------------------------------------Epidemic size---------------------------------------------------*/

void Epid_siz(Matdoub Coo,Matdoub Selke,Matdoub Param, int N, double Tmax,  double dt, Vecdoub testim, int nr,int dtype,int nt,Vecdoub map,Vecdoub L,Vecdoub l,Vecint &M,int Nepidem){
  M.resize(Nepidem); Vecdoub selke(N);
  for(int i=0;i<Nepidem;i++){
    M[i]=SI_mimulation_test2(Coo,Selke,Param,N,Tmax,dt,testim,nr,dtype,nt, map,L,l,i);
  }
  
}


