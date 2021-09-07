#include <stdio.h>
#include <assert.h>
#include "bspline.h"
#include "common.h"

/*
  reference:
    A smooth particle mesh Ewald method, J. Chem. Phys. 103, 8577 (1995);

  bspline_1d and bspline_2d are kept in simple form.
  bspline_3d has some extra optimization and no other extra for the higher-dimention ones
*/

__global__ void bspline_1d(double *data, double min_x1, double max_x1,
  int num_x1, double x1, int order, double* o_value, 
  double* o_d1, bool f_d1, double* o_d2, bool f_d2)
/*
  num_x1: number of data points
  min_x1: x-coordinate of the first data point
  max_x1: x-coordinate of the last data point
  order: B-spline order
  o_value: fitting value
  o_d1: first-derivative value
  f_d1: flag for calculating first-derivative value
  o_d2: second-derivative value
  f_d2: flag for calculating second-derivative value
*/
{
  int numbin_x1 = num_x1 - 1; // numbin_x1: number of bins, i.e. number of data points minus 1 (data points are on the grids)
  double x1scaled = (x1 - min_x1) / (max_x1 - min_x1) * numbin_x1;

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1) {
    printf("Error: Out of range (bspline_1d)");
    __trap();
  }

  // If index is nearby the boundary, adjust accordingly.
  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  /*
    Mn(order, ...) has domain [0, order] with center (peak) at order/2.0,
    so later on the calculated relative distance to the grid point (x1scaled - i)
    need to be shifted to correspond to the peak.
  */
  x1scaled += order/2.0; 

  for(int i = x1grid0; i <= x1grid1; i++){
    double datapoint = *(data + i);

    value += datapoint * Mn(order, x1scaled - i);

    if(f_d1){
      d1 += datapoint * dMn(order, x1scaled - i) 
        * (numbin_x1 / (max_x1 - min_x1));
    }

    if(f_d2){
      d2 += datapoint * dMn2(order, x1scaled - i)
        * (numbin_x1 / (max_x1 - min_x1)) * (numbin_x1 / (max_x1 - min_x1));
    }
  }

  *o_value = value;
  *o_d1 = d1;
  *o_d2 = d2;
}

__global__ void bspline_2d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  int order, double* o_value, double* o_d1, double* o_d2,
  bool f_d)
/*
  o_value: fitting value
  o_d1: partial first-derivative value along x1 direction
  o_d2: partial first-derivative value along x2 direction (not to be confused with the one in bspline_1d)
  f_d1: flag for calculating partial first-derivative value
*/    
{
  int numbin_x1 = num_x1 - 1;
  int numbin_x2 = num_x2 - 1;

  double x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1;
  double x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2;

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  int x2grid0 = int(x2scaled - order/2.0) + 1;
  int x2grid1 = int(x2scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1 || x2grid1 < 0 || x2grid0 > numbin_x2) {
    printf("Error: Out of range (bspline_2d)");
    __trap();
  }

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 

  for(int i = x1grid0; i <= x1grid1; i++){
    for(int j = x2grid0; j <= x2grid1; j++){
      double datapoint = *(data+i*num_x2+j);

      value += datapoint*Mn(order, x1scaled - i)*Mn(order, x2scaled - j);

      if(f_d){
        d1 += datapoint*dMn(order, x1scaled - i)*Mn(order, x2scaled - j) 
          *(numbin_x1/(max_x1-min_x1));

        d2 += datapoint*Mn(order, x1scaled - i)*dMn(order, x2scaled - j) 
          *(numbin_x2/(max_x2-min_x2));
      }
    }
  }

  *o_value = value;
  *o_d1 = d1;
  *o_d2 = d2;
}
  
__global__ void bspline_3d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double* o_value, double* o_d1, double* o_d2, double* o_d3,
  bool f_d)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  int numbin_x1 = num_x1 - 1;
  int numbin_x2 = num_x2 - 1;
  int numbin_x3 = num_x3 - 1;
  
  double x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1;
  double x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2;
  double x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  int x2grid0 = int(x2scaled - order/2.0) + 1;
  int x2grid1 = int(x2scaled + order/2.0);

  int x3grid0 = int(x3scaled - order/2.0) + 1;
  int x3grid1 = int(x3scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1 
    || x2grid1 < 0 || x2grid0 > numbin_x2
    || x3grid1 < 0 || x3grid0 > numbin_x3) {
    if(idx == 0) printf("Error: Out of range (bspline_3d)");
    __trap();
  }

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  if(x3grid1 > numbin_x3) x3grid1 = numbin_x3;
  if(x3grid0 < 0) x3grid0 = 0;

  int r1 = x1grid1 - x1grid0 + 1;
  int r2 = x2grid1 - x2grid0 + 1;
  int r3 = x3grid1 - x3grid0 + 1;

  int i = idx/(r2*r3);
  int j = (idx-i*r2*r3)/r3;
  int k = idx-i*r2*r3-j*r3;

  i += x1grid0;
  j += x2grid0;
  k += x3grid0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;

  /* only in-bound threads get datapoint 
     out-bound threads remain active for shuffle but feed value 0
  */
  double datapoint = 0.0;
  if(idx < r1*r2*r3)   
    datapoint = *(data+i*num_x2*num_x3
    +j*num_x3+k);

  /*
    breaking them down here since they are reusable for the first-derivative if requested
  */ 
  double Mn1a = Mn(order-1, x1scaled - i);
  double Mn1b = Mn(order-1, x1scaled - i - 1);
  double Mn2a = Mn(order-1, x2scaled - j);
  double Mn2b = Mn(order-1, x2scaled - j - 1);
  double Mn3a = Mn(order-1, x3scaled - k);
  double Mn3b = Mn(order-1, x3scaled - k - 1);

  double Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b;
  double Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b;
  double Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b;

  double value = datapoint*Mn1*Mn2*Mn3;

  value = blockReduceSum(value);

  if (threadIdx.x==0){
    *(o_value+blockIdx.x) = value;
  }

  if(f_d){
    double dMn1 = Mn1a - Mn1b;
    double dMn2 = Mn2a - Mn2b;
    double dMn3 = Mn3a - Mn3b;

    double d1 = datapoint*dMn1*Mn2*Mn3
      *(numbin_x1/(max_x1-min_x1));

    double d2 = datapoint*Mn1*dMn2*Mn3
      *(numbin_x2/(max_x2-min_x2));

    double d3 = datapoint*Mn1*Mn2*dMn3
      *(numbin_x3/(max_x3-min_x3));

    d1 = blockReduceSum(d1);
    d2 = blockReduceSum(d2);
    d3 = blockReduceSum(d3);

    if (threadIdx.x==0){
      *(o_d1+blockIdx.x) = d1;
      *(o_d2+blockIdx.x) = d2;
      *(o_d3+blockIdx.x) = d3;
    }
  }
}

__global__ void bspline_4d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double* o_value, double* o_d1, double* o_d2, double* o_d3, double* o_d4,
  bool f_d)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;

  int numbin_x1 = num_x1 - 1;
  int numbin_x2 = num_x2 - 1;
  int numbin_x3 = num_x3 - 1;
  int numbin_x4 = num_x4 - 1;
    
  double x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1;
  double x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2;
  double x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3;
  double x4scaled = (x4-min_x4)/(max_x4-min_x4)*numbin_x4;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  int x2grid0 = int(x2scaled - order/2.0) + 1;
  int x2grid1 = int(x2scaled + order/2.0);

  int x3grid0 = int(x3scaled - order/2.0) + 1;
  int x3grid1 = int(x3scaled + order/2.0);

  int x4grid0 = int(x4scaled - order/2.0) + 1;
  int x4grid1 = int(x4scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1 
    || x2grid1 < 0 || x2grid0 > numbin_x2
    || x3grid1 < 0 || x3grid0 > numbin_x3
    || x4grid1 < 0 || x4grid0 > numbin_x4) {
    if(idx == 0) printf("Error: Out of range (bspline_4d)");
    __trap();
  }

  /*
  assert(x1grid0 >= 0 && x1grid1 <= numbin_x1 
    && x2grid0 >= 0 && x2grid1 <= numbin_x2
    && x3grid0 >= 0 && x3grid1 <= numbin_x3
    && x4grid0 >= 0 && x4grid1 <= numbin_x4); 
  */

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  if(x3grid1 > numbin_x3) x3grid1 = numbin_x3;
  if(x3grid0 < 0) x3grid0 = 0;

  if(x4grid1 > numbin_x4) x4grid1 = numbin_x4;
  if(x4grid0 < 0) x4grid0 = 0;

  int r1 = x1grid1 - x1grid0 + 1;
  int r2 = x2grid1 - x2grid0 + 1;
  int r3 = x3grid1 - x3grid0 + 1;
  int r4 = x4grid1 - x4grid0 + 1;

  int i = idx/(r2*r3*r4);
  int j = (idx-i*r2*r3*r4)/(r3*r4);
  int k = (idx-i*r2*r3*r4-j*r3*r4)/r4;
  int m = idx-i*r2*r3*r4-j*r3*r4-k*r4;

  i += x1grid0;
  j += x2grid0;
  k += x3grid0;
  m += x4grid0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;
  x4scaled += order/2.0;

  /* only in-bound threads get datapoint 
     out-bound threads remain active for shuffle but feed value 0
  */
  double datapoint = 0.0;
  if(idx < r1*r2*r3*r4) 
    datapoint = *(data+i*num_x2*num_x3*num_x4
    +j*num_x3*num_x4+k*num_x4+m);

  double Mn1a = Mn(order-1, x1scaled - i);
  double Mn1b = Mn(order-1, x1scaled - i - 1);
  double Mn2a = Mn(order-1, x2scaled - j);
  double Mn2b = Mn(order-1, x2scaled - j - 1);
  double Mn3a = Mn(order-1, x3scaled - k);
  double Mn3b = Mn(order-1, x3scaled - k - 1);
  double Mn4a = Mn(order-1, x4scaled - m);
  double Mn4b = Mn(order-1, x4scaled - m - 1);

  double Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b;
  double Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b;
  double Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b;
  double Mn4 = (x4scaled - m)/(order-1)*Mn4a + (order-(x4scaled - m))/(order-1)*Mn4b;

  double value = datapoint*Mn1*Mn2*Mn3*Mn4;

  value = blockReduceSum(value);

  if (threadIdx.x==0){
    *(o_value+blockIdx.x) = value;
  }

  if(f_d){
    double dMn1 = Mn1a - Mn1b;
    double dMn2 = Mn2a - Mn2b;
    double dMn3 = Mn3a - Mn3b;
    double dMn4 = Mn4a - Mn4b;

    double d1 = datapoint*dMn1*Mn2*Mn3*Mn4
      *(numbin_x1/(max_x1-min_x1));

    double d2 = datapoint*Mn1*dMn2*Mn3*Mn4
      *(numbin_x2/(max_x2-min_x2));

    double d3 = datapoint*Mn1*Mn2*dMn3*Mn4
      *(numbin_x3/(max_x3-min_x3));

    double d4 = datapoint*Mn1*Mn2*Mn3*dMn4 
      *(numbin_x4/(max_x4-min_x4));

    d1 = blockReduceSum(d1);
    d2 = blockReduceSum(d2);
    d3 = blockReduceSum(d3);
    d4 = blockReduceSum(d4);

    if (threadIdx.x==0){
      *(o_d1+blockIdx.x) = d1;
      *(o_d2+blockIdx.x) = d2;
      *(o_d3+blockIdx.x) = d3;
      *(o_d4+blockIdx.x) = d4;
    }
  }

}

__global__ void bspline_5d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double* o_value, double* o_d1, double* o_d2, double* o_d3, double* o_d4, double* o_d5,
  bool f_d)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
 
  int numbin_x1 = num_x1 - 1;
  int numbin_x2 = num_x2 - 1;
  int numbin_x3 = num_x3 - 1;
  int numbin_x4 = num_x4 - 1;
  int numbin_x5 = num_x5 - 1;

  double x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1;
  double x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2;
  double x3scaled = (x3-min_x3)/(max_x3-min_x3)*numbin_x3;
  double x4scaled = (x4-min_x4)/(max_x4-min_x4)*numbin_x4;
  double x5scaled = (x5-min_x5)/(max_x5-min_x5)*numbin_x5;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  int x2grid0 = int(x2scaled - order/2.0) + 1;
  int x2grid1 = int(x2scaled + order/2.0);

  int x3grid0 = int(x3scaled - order/2.0) + 1;
  int x3grid1 = int(x3scaled + order/2.0);

  int x4grid0 = int(x4scaled - order/2.0) + 1;
  int x4grid1 = int(x4scaled + order/2.0);

  int x5grid0 = int(x5scaled - order/2.0) + 1;
  int x5grid1 = int(x5scaled + order/2.0);


  if(x1grid1 < 0 || x1grid0 > numbin_x1 
    || x2grid1 < 0 || x2grid0 > numbin_x2
    || x3grid1 < 0 || x3grid0 > numbin_x3
    || x4grid1 < 0 || x4grid0 > numbin_x4
    || x5grid1 < 0 || x5grid0 > numbin_x5) {
    if(idx == 0) printf("Error: Out of range (bspline_5d)");
    __trap();
  }

  /*
  assert(x1grid0 >= 0 && x1grid1 <= numbin_x1 
    && x2grid0 >= 0 && x2grid1 <= numbin_x2
    && x3grid0 >= 0 && x3grid1 <= numbin_x3
    && x4grid0 >= 0 && x4grid1 <= numbin_x4
    && x5grid0 >= 0 && x5grid1 <= numbin_x5); 
  */

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  if(x3grid1 > numbin_x3) x3grid1 = numbin_x3;
  if(x3grid0 < 0) x3grid0 = 0;

  if(x4grid1 > numbin_x4) x4grid1 = numbin_x4;
  if(x4grid0 < 0) x4grid0 = 0;

  if(x5grid1 > numbin_x5) x5grid1 = numbin_x5;
  if(x5grid0 < 0) x5grid0 = 0;

  int r1 = x1grid1 - x1grid0 + 1;
  int r2 = x2grid1 - x2grid0 + 1;
  int r3 = x3grid1 - x3grid0 + 1;
  int r4 = x4grid1 - x4grid0 + 1;
  int r5 = x5grid1 - x5grid0 + 1;

  int i = idx/(r2*r3*r4*r5);
  int j = (idx-i*r2*r3*r4*r5)/(r3*r4*r5);
  int k = (idx-i*r2*r3*r4*r5-j*r3*r4*r5)/(r4*r5);
  int m = (idx-i*r2*r3*r4*r5-j*r3*r4*r5-k*r4*r5)/r5;
  int n = idx-i*r2*r3*r4*r5-j*r3*r4*r5-k*r4*r5-m*r5;

  i += x1grid0;
  j += x2grid0;
  k += x3grid0;
  m += x4grid0;
  n += x5grid0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;
  x4scaled += order/2.0;
  x5scaled += order/2.0;

  /* only in-bound threads get datapoint 
     out-bound threads remain active for shuffle but feed value 0
  */
  double datapoint = 0.0;
  if(idx < r1*r2*r3*r4*r5) 
    datapoint = *(data+i*num_x2*num_x3*num_x4*num_x5
    +j*num_x3*num_x4*num_x5+k*num_x4*num_x5+m*num_x5+n);

  double Mn1a = Mn(order-1, x1scaled - i);
  double Mn1b = Mn(order-1, x1scaled - i - 1);
  double Mn2a = Mn(order-1, x2scaled - j);
  double Mn2b = Mn(order-1, x2scaled - j - 1);
  double Mn3a = Mn(order-1, x3scaled - k);
  double Mn3b = Mn(order-1, x3scaled - k - 1);
  double Mn4a = Mn(order-1, x4scaled - m);
  double Mn4b = Mn(order-1, x4scaled - m - 1);
  double Mn5a = Mn(order-1, x5scaled - n);
  double Mn5b = Mn(order-1, x5scaled - n - 1);


  double Mn1 = (x1scaled - i)/(order-1)*Mn1a + (order-(x1scaled - i))/(order-1)*Mn1b;
  double Mn2 = (x2scaled - j)/(order-1)*Mn2a + (order-(x2scaled - j))/(order-1)*Mn2b;
  double Mn3 = (x3scaled - k)/(order-1)*Mn3a + (order-(x3scaled - k))/(order-1)*Mn3b;
  double Mn4 = (x4scaled - m)/(order-1)*Mn4a + (order-(x4scaled - m))/(order-1)*Mn4b;
  double Mn5 = (x5scaled - n)/(order-1)*Mn5a + (order-(x5scaled - n))/(order-1)*Mn5b;

  double value = datapoint*Mn1*Mn2*Mn3*Mn4*Mn5;

  value = blockReduceSum(value);

  if (threadIdx.x==0){
    *(o_value+blockIdx.x) = value;
  }  

  if(f_d){
    double dMn1 = Mn1a - Mn1b;
    double dMn2 = Mn2a - Mn2b;
    double dMn3 = Mn3a - Mn3b;
    double dMn4 = Mn4a - Mn4b;
    double dMn5 = Mn5a - Mn5b;

    double d1 = datapoint*dMn1*Mn2*Mn3*Mn4*Mn5
      *(numbin_x1/(max_x1-min_x1));

    double d2 = datapoint*Mn1*dMn2*Mn3*Mn4*Mn5
      *(numbin_x2/(max_x2-min_x2));

    double d3 = datapoint*Mn1*Mn2*dMn3*Mn4*Mn5
      *(numbin_x3/(max_x3-min_x3));

    double d4 = datapoint*Mn1*Mn2*Mn3*dMn4*Mn5 
      *(numbin_x4/(max_x4-min_x4));

    double d5 = datapoint*Mn1*Mn2*Mn3*Mn4*dMn5 
      *(numbin_x5/(max_x5-min_x5));

    d1 = blockReduceSum(d1);
    d2 = blockReduceSum(d2);
    d3 = blockReduceSum(d3);
    d4 = blockReduceSum(d4);
    d5 = blockReduceSum(d5);

    if (threadIdx.x==0){
      *(o_d1+blockIdx.x) = d1;
      *(o_d2+blockIdx.x) = d2;
      *(o_d3+blockIdx.x) = d3;
      *(o_d4+blockIdx.x) = d4;
      *(o_d5+blockIdx.x) = d5;
    }  
  }

}

void bspline_5d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4, double& r_d5,
  bool f_d)
{
  int size = order*order*order*order*order; 
  int n_block = (size+BLOCK_SIZE-1) / BLOCK_SIZE;

  double *value,*d1,*d2,*d3,*d4,*d5;
  CHECK(cudaMallocManaged((void **)&value, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d1, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d2, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d3, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d4, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d5, n_block*sizeof(double)));

  bspline_5d<<< n_block, BLOCK_SIZE >>>(data, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    min_x5, max_x5, num_x5, x5,
    order, value, d1, d2, d3, d4, d5, f_d);
  cudaDeviceSynchronize();

  r_value = 0.0;
  r_d1 = 0.0;
  r_d2 = 0.0;
  r_d3 = 0.0;
  r_d4 = 0.0;
  r_d5 = 0.0;
  for(int i = 0; i < n_block; i++){
    r_value += value[i];
    r_d1 += d1[i];
    r_d2 += d2[i];
    r_d3 += d3[i];
    r_d4 += d4[i];
    r_d5 += d5[i];
  }

  CHECK(cudaFree(value));
  CHECK(cudaFree(d1));
  CHECK(cudaFree(d2));
  CHECK(cudaFree(d3));
  CHECK(cudaFree(d4));
  CHECK(cudaFree(d5));
}

void bspline_4d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4,
  bool f_d)
{
  int size = order*order*order*order; 
  int n_block = (size+BLOCK_SIZE-1) / BLOCK_SIZE;

  double *value,*d1,*d2,*d3,*d4;
  CHECK(cudaMallocManaged((void **)&value, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d1, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d2, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d3, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d4, n_block*sizeof(double)));

  bspline_4d<<< n_block, BLOCK_SIZE >>>(data, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    order, value, d1, d2, d3, d4, f_d);
  cudaDeviceSynchronize();

  r_value = 0.0;
  r_d1 = 0.0;
  r_d2 = 0.0;
  r_d3 = 0.0;
  r_d4 = 0.0;
  for(int i = 0; i < n_block; i++){
    r_value += value[i];
    r_d1 += d1[i];
    r_d2 += d2[i];
    r_d3 += d3[i];
    r_d4 += d4[i];
  }

  CHECK(cudaFree(value));
  CHECK(cudaFree(d1));
  CHECK(cudaFree(d2));
  CHECK(cudaFree(d3));
  CHECK(cudaFree(d4));
}

void bspline_3d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3,
  bool f_d)
{
  int size = order*order*order*order; 
  int n_block = (size+BLOCK_SIZE-1) / BLOCK_SIZE;

  double *value,*d1,*d2,*d3;
  CHECK(cudaMallocManaged((void **)&value, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d1, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d2, n_block*sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d3, n_block*sizeof(double)));

  bspline_3d<<< n_block, BLOCK_SIZE >>>(data, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    order, value, d1, d2, d3, f_d);
  cudaDeviceSynchronize();

  r_value = 0.0;
  r_d1 = 0.0;
  r_d2 = 0.0;
  r_d3 = 0.0;
  for(int i = 0; i < n_block; i++){
    r_value += value[i];
    r_d1 += d1[i];
    r_d2 += d2[i];
    r_d3 += d3[i];
  }

  CHECK(cudaFree(value));
  CHECK(cudaFree(d1));
  CHECK(cudaFree(d2));
  CHECK(cudaFree(d3));
}
