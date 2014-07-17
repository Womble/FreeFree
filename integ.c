//compile me with gcc -fPIC -shared -O3 -o integ.so integ.c
//for use as a shared library within python
#include "math.h"
#include "stdio.h"

double integrate (double *column, double cellLength, int size){
  double intensity=0, expDt;
  int i,j;
  for (i=0;i<size/3;i++){
    /* if (fabs(column[2*i+1])>1e-50 && fabs(column[2*i])>1e-50){ */
    /*   fprintf(stdout, "%.3e %.3e %.3e\n",intensity,column[2*i+1],column[2*i]); */
    /* } */
    for (j=0;j<1000;j++){
      intensity=intensity+cellLength/1000*(-column[2*i+1]*intensity+column[2*i+1]);
    }
  }
  return intensity;
}
