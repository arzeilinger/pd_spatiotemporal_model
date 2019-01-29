 
//Indx of an elment
template <class T,class U>
Int which(T &rhs, U key){
	for (Int i=0;i<rhs.size();i++) {
	 if(rhs[i]==key){
	   return i;
	 }
	}
}

// Check whether an element is in a vector
template <class T,class U>
Bool Check(T &a, U &b){
bool B=true;
int s=0;
 for(int i=0;i<a.size();i++){
  if(a[i]==b){
   s+=1;
   break;
  }
 }
 if(s==0){
  return false;
 }
 else{
  return B;
 }
}
// Number of element two vectors have in common
template <class T>
Int Num_element(T &a, T &b){
 Int s=0;
  for(int i=0; i<a.size();i++){
   for(int j=0; j<b.size();j++){
    if(a[i]==b[j]){
     s+=1;
     break;
    }
  }
  }
  return s;
}
//Delete an integer b from a vector of integer a ie a[!b]
template <class T,class U>
void Delet(T &a,U &b,T &c){
 int k= which(a,b);
 for(int i=0;i<a.size();i++){
  if(i<k){
   c[i]=a[i];
  }
  if(i>k){
   c[i-1]=a[i];
  }
 }
}

//delete b elments from a vector ie a[!b]
template <class T>
void Delet_vec(T &a,T &b,T &c){
int k=0;
 for(int j=0;j<a.size();j++){
  for(int i=0;i<b.size();i++){
  if(a[i]==b[j]){
   c[k]=a[i];
   k+=1;
   break;
  }
 }
 }
}

// Get a subset of a matrix by row
void mat_subset(MatDoub &a, Int &n,VecInt &c, MatDoub &b){
 for(int i=0;i<n;i++){
  for(int j=0;j<c.size();j++){
   b[i][c[j]]=a[i][c[j]];
  }
 }
}

//Distance
Doub dist(VecDoub &a, MatDoub &b){
 Doub s=0;
 for(int i=0; i< b.nrows();i++){
  s+=sqrt(pow(a[0]-b[i][0],2)+ pow(a[1]-b[i][1],2));

 }
return s;
}

Doub dist1(VecDoub &a, VecDoub &b){
 Doub s=sqrt(pow(a[0]-b[0],2)+ pow(a[0]-b[1],2));
return s;
}

// Distance between a point and a set of point
void dist2(VecDoub &a, MatDoub &b,VecDoub &c){
for(int i=0;i<b.nrows();i++){
  VecDoub e(2);
  e[0]=b[i][0];
  e[1]=b[i][1];
  c[i]=dist1(a,e);
}

}

// copy vector
void copy_vec(VecDoub &a,VecDoub &b){
Int n=a.size();
for(int i=0;i<n;i++){
 b[i]=a[i];
}
}
