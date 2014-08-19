//compile me with gcc -fPIC -shared -O3 -o integ.so integ.c
//for use as a shared library within python
#include "math.h"
#include "stdio.h"

double integrate (double *column, double cellLength, int size){
  //column is a 1D array for the form [dt_0,S_0,dt_1,S_1,dt_2...,dt_size/2,S_size/2]
  double intensity=0, dt,S;
  int i,j,n;
  for (i=0;i<size/2;i++){
    if (column[2*i]>0){ //we're not going to even attempt negative optical depths, assuming they are erronious rather than masing
      if (column[2*i]>10) intensity=column[2*i+1]; //for cells of high optical depth set intensity to saturated value
      else{
      	n=(20*column[2*i]);
      	if (n<1)   n=1;
      	dt=column[2*i]/n;
      	for (j=0;j<n;j++){ //#integrate in steps of tau~=1/20 i.e intensity drop of ~5%
      	  //	printf("%d: %.3e -> ",i,intensity);
      	  intensity+=dt*(column[2*i+1]-intensity);
      	  if(intensity<0) intensity=0; //this is probably unnecessary, but im being paranoid atm
      	  //printf("%.3e %.3e -> %.3e\n",dt,column[2*i+1], intensity);
      	}
      }
      //intensity=intensity*exp(-column[2*i])+column[2*i+1]*column[2*i];
    }
  }
  return intensity;
}
