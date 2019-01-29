#include <fstream>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
using namespace std;

typedef vector<int> Vecint;
typedef vector<double> Vecdoub;
typedef vector<vector<int> > Matint;
typedef vector<vector<double> > Matdoub;

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

template<class T>
inline T SQR(const T a) {return a*a;}

template<class T>
inline const T &MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}
        
        template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}


//template <typename T, size_t size>
//T* begin(T (&c)[size])
//{
//    return &c[0];
//}

//template <typename T, size_t size>
//T* end(T (&c)[size])
//{
//    return &c[0] + size;
//}


//Indx of an elment
template <class T,class U>
inline int which(T rhs, U key){
	for (int i=0;i<rhs.size();i++) {
	 if(rhs[i]==key){
	   return i;
	   break;
	 }
	}
}

// Check whether an element is in a vector
template <class T,class U >
inline bool Check(T a, U b){
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
int Num_element(T &a, T &b){
 int s=0;
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
inline void Delet(T a,U b,T &c){
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
inline void Delet_vec(T &a,T b){
int k=0;
 //for(int j=0;j<a.size();j++){
  for(int i=0;i<b.size();i++){
  //if(a[i]==b[j]){
   a.erase(remove(a.begin(),a.end(),b[i]),a.end());
 //}
 }
 //}
}
//Get a subset of a vector
template <class T,class U>
inline void Vec_subset(T a,U &b,int n,int m){
if(m<=n){
cout<<"m must be greate than n"<<endl;
}
else{
 b.resize(m-n);

  for(int i=n;i<m;i++){
   b[i]=a[i];
  }
}

}
// Get a subset of a matrix by row
inline void mat_subset(Matdoub a, int n,Vecint c, Matdoub &b){
b.resize(n);
int m=c.size();

for(int i=0;i<n;i++){
 b[i].resize(m);
}
 for(int i=0;i<n;i++){
  for(int j=0;j<m;j++){
   b[i][j]=a[i][c[j]];
//   cout<<b[i][j]<<endl;
  }
 }
}

// Seletion of some rows of a matrix

inline void Selec(Matdoub a, Vecint v,Vecint c, Matdoub &b){
int n=v.size();

b.resize(n);
int m=c.size();

for(int i=0;i<n;i++){
 b[i].resize(m);
}
 for(int i=0;i<n;i++){
  for(int j=0;j<m;j++){
   b[i][j]=a[v[i]-1][c[j]];
//   cout<<b[i][j]<<endl;
  }
 }
}

//Distance
template <class T,class U>
double dist(T &a, U &b){
int nrows= sizeof(b)/sizeof(b[0]);
 double s=0;
 for(int i=0; i< nrows;i++){
  s+=sqrt(pow(a[0]-b[i][0],2)+ pow(a[1]-b[i][1],2));

 }
return s;
}
template <class T>
double dist1(T a, T b){
 double s=sqrt(pow(a[0]-b[0],2)+ pow(a[1]-b[1],2));
return s;
}

// Distance between a point and a set of point
template <class T,class U>
void dist2(T a, U b,T &c,int n){
//n is the size of the vectors
//cout<<nrows<<endl;
for(int i=0;i<n;i++){
  T e(2);
  e[0]=b[i][0];
  e[1]=b[i][1];
  c[i]=dist1(a,e);
}

}
// Function to pass a file and return it into an array
template <class T,class U>
void file_mat(T &infile,int nrow,int ncol, U &arr,string filename)
{
//  string  filename;
//  cin >>filename;
  infile.open(filename.c_str());
  arr.resize(nrow);
  for(int i=0;i<nrow;i++){
   arr[i].resize(ncol); 
  }
 for (int i = 0; i < nrow; i++){
       for (int j = 0; j < ncol; j++){
      infile>> arr[i][j] ;
    }
            
    }
       infile.close(); 
}

//Function to pass file and return a vector
template <class T,class U>
void file_vec(T &infile,int nrow,U &arr,string filename)
{
//  string  filename;
//  cin >>filename;
  infile.open(filename.c_str());
  arr.resize(nrow);
 for (int i = 0; i < nrow; i++){
      infile>> arr[i] ;
    }
       infile.close(); 
}


// kernel function
double kernel(double d,double alpha,double dtype ){
const double pi = 3.1415926535897;
double res=1/(pi*(2*alpha+dtype));
if(d<=dtype){
 res*=1/dtype;
}
else{
  res*=1/d*exp(-(d-dtype)/alpha);
}
return res;
}

// sum of all kernels
double sum_kernel(Vecdoub d,double alpha,int dtype ){
double res=0;
for(int i=0;i<d.size();i++){
 //for(int j=0; j<d.ncols();j++){
  res+=kernel(d[i],alpha,dtype);
//}
}
return res;
}
//function to write a data
template <class T,class U>
void ou_file(T &ofile,int nrow,int ncol, U &arr,string filename)
{
  ofile.open(filename);
    arr.resize(nrow);
    for(int i=0;i<nrow;i++){
     arr[i].resize(ncol); 
    }
    for (int i = 0; i < nrow; i++){
       for (int j = 0; j < ncol; j++){
      ofile<< arr[i][j]<<" " ;
    }
    ofile<<endl;        
    }
       ofile.close(); 
}
template <class T,class U>
void ou_file_vec(T &ofile,int ncol, U &arr,string filename)
{
  ofile.open(filename);
    arr.resize(ncol);
    for(int i=0;i<ncol;i++){
     ofile<<arr[i]<<" ";
    }
     ofile<<endl;  
        ofile.close(); 
}

// Data form
void data(Matdoub Coo, Vecdoub Inf_time, Vecint Indx_ind, Matdoub &M ,int NI){
//int n=Inf_time.size();
M.resize(NI);
for(int i=0;i<NI;i++){
  M[i].resize(4);
}
int i,k;
Vecdoub t_sort(Inf_time);
sort(t_sort.begin(),t_sort.end());
//for(int i=0;i<200;i++){
//  cout<<t_sort[i]<<" ";
//}
//cout<<endl;
//cout<<t_sort[46]<<endl;
for(i=0;i<NI;i++){ 
 k=which(Inf_time,t_sort[i]);
// cout<<t_sort[i]<<endl;
 M[i][0]=t_sort[i];
 M[i][1]=Coo[k][0];
 M[i][2]=Coo[k][1];
 M[i][3]=k;
}

}

