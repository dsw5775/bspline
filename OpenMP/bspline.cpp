#include <stdio.h>
#include <stdlib.h>
#include "bspline.h"

/*
  reference:
    A smooth particle mesh Ewald method, J. Chem. Phys. 103, 8577 (1995);

  bspline_1d and bspline_2d are kept in simple form.
  bspline_3d has some extra optimization and no other extra for the higher-dimention ones
*/
 
void bspline_1d(double *data, double min_x1, double max_x1,
  int num_x1, double x1, int order, double& r_value, 
  double& r_d1, bool f_d1, double& r_d2, bool f_d2)
/*
  num_x1: number of data points
  min_x1: x-coordinate of the first data point
  max_x1: x-coordinate of the last data point
  order: B-spline order
  r_value: fitting value
  r_d1: first-derivative value
  f_d1: flag for calculating first-derivative value
  r_d2: second-derivative value
  f_d2: flag for calculating second-derivative value
*/
{
  int numbin_x1 = num_x1 - 1; // numbin_x1: number of bins, i.e. number of data points minus 1 (data points are on the grids)
  double x1scaled = (x1 - min_x1) / (max_x1 - min_x1) * numbin_x1;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1) {
    printf("out of range (bspline_1d)\n"); 
    exit(-1);
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

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;

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

  r_value = value;
  r_d1 = d1;
  r_d2 = d2;
}

void bspline_2d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  int order, double& r_value, double& r_d1, double& r_d2, 
  bool f_d)
/*
  r_value: fitting value
  r_d1: partial first-derivative value along x1 direction
  r_d2: partial first-derivative value along x2 direction (not to be confused with the one in bspline_1d)
  f_d1: flag for calculating partial first-derivative value
*/  
{
  int numbin_x1 = num_x1 - 1;
  int numbin_x2 = num_x2 - 1;

  double x1scaled = (x1-min_x1)/(max_x1-min_x1)*numbin_x1;
  double x2scaled = (x2-min_x2)/(max_x2-min_x2)*numbin_x2;

  int x1grid0 = int(x1scaled - order/2.0) + 1;
  int x1grid1 = int(x1scaled + order/2.0);

  int x2grid0 = int(x2scaled - order/2.0) + 1;
  int x2grid1 = int(x2scaled + order/2.0);

  if(x1grid1 < 0 || x1grid0 > numbin_x1 || x2grid1 < 0 || x2grid0 > numbin_x2) {
    printf("out of range (bspline_2d)\n"); 
    exit(-1);
  }

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;

  for(int i = x1grid0; i <= x1grid1; i++){
    for(int j = x2grid0; j <= x2grid1; j++){
      double datapoint = *(data+i*num_x2+j);

      double Mn1 = Mn(order, x1scaled - i);
      double Mn2 = Mn(order, x2scaled - j);

      value += datapoint*Mn1*Mn2;

      if(f_d){
        double dMn1 = dMn(order, x1scaled - i);
        double dMn2 = dMn(order, x2scaled - j);

        d1 += datapoint*dMn1*Mn2 
          *(numbin_x1/(max_x1-min_x1));

        d2 += datapoint*Mn1*dMn2 
          *(numbin_x2/(max_x2-min_x2));
      }
    }
  }

  r_value = value;
  r_d1 = d1;
  r_d2 = d2;
}
  
void bspline_3d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3,
  bool f_d)
{
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
    printf("out of range (bspline_3d)\n"); 
    exit(-1);
  }

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  if(x3grid1 > numbin_x3) x3grid1 = numbin_x3;
  if(x3grid0 < 0) x3grid0 = 0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  double d3 = 0.0;

  #pragma omp parallel for collapse(3) reduction(+:value,d1,d2,d3)
  for(int i = x1grid0; i <= x1grid1; i++){
    for(int j = x2grid0; j <= x2grid1; j++){
      for(int k = x3grid0; k <= x3grid1; k++){
        double datapoint = *(data+i*num_x2*num_x3+j*num_x3+k);

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
 
        value += datapoint*Mn1*Mn2*Mn3;
  
        if(f_d){
          double dMn1 = Mn1a - Mn1b;
          double dMn2 = Mn2a - Mn2b;
          double dMn3 = Mn3a - Mn3b;

          d1 += datapoint*dMn1*Mn2*Mn3 
            *(numbin_x1/(max_x1-min_x1));
  
          d2 += datapoint*Mn1*dMn2*Mn3
            *(numbin_x2/(max_x2-min_x2));

          d3 += datapoint*Mn1*Mn2*dMn3 
            *(numbin_x3/(max_x3-min_x3));
        }
      }
    }
  }

  r_value = value;
  r_d1 = d1;
  r_d2 = d2;
  r_d3 = d3;
}

void bspline_4d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4,
  bool f_d)
{
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
    printf("out of range (bspline_4d)\n"); 
    exit(-1);
  }

  if(x1grid1 > numbin_x1) x1grid1 = numbin_x1;
  if(x1grid0 < 0) x1grid0 = 0;

  if(x2grid1 > numbin_x2) x2grid1 = numbin_x2;
  if(x2grid0 < 0) x2grid0 = 0;

  if(x3grid1 > numbin_x3) x3grid1 = numbin_x3;
  if(x3grid0 < 0) x3grid0 = 0;

  if(x4grid1 > numbin_x4) x4grid1 = numbin_x4;
  if(x4grid0 < 0) x4grid0 = 0;

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;
  x4scaled += order/2.0;

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  double d3 = 0.0;
  double d4 = 0.0;

  #pragma omp parallel for collapse(4) reduction(+:value,d1,d2,d3,d4) 
  for(int i = x1grid0; i <= x1grid1; i++){
    for(int j = x2grid0; j <= x2grid1; j++){
      for(int k = x3grid0; k <= x3grid1; k++){
        for(int m = x4grid0; m <= x4grid1; m++){
          double datapoint = *(data+i*num_x2*num_x3*num_x4+j*num_x3*num_x4+k*num_x4+m);

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

          value += datapoint*Mn1*Mn2*Mn3*Mn4;
    
          if(f_d){
            double dMn1 = Mn1a - Mn1b;
            double dMn2 = Mn2a - Mn2b;
            double dMn3 = Mn3a - Mn3b;
            double dMn4 = Mn4a - Mn4b;

            d1 += datapoint*dMn1*Mn2*Mn3*Mn4
              *(numbin_x1/(max_x1-min_x1));
    
            d2 += datapoint*Mn1*dMn2*Mn3*Mn4
              *(numbin_x2/(max_x2-min_x2));
  
            d3 += datapoint*Mn1*Mn2*dMn3*Mn4 
              *(numbin_x3/(max_x3-min_x3));

            d4 += datapoint*Mn1*Mn2*Mn3*dMn4 
              *(numbin_x4/(max_x4-min_x4));
          }
        }
      }
    }
  }

  r_value = value;
  r_d1 = d1;
  r_d2 = d2;
  r_d3 = d3;
  r_d4 = d4;
}

void bspline_5d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4, double& r_d5,
  bool f_d)
{
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
    printf("out of range (bspline_5d)\n"); 
    exit(-1);
  }

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

  x1scaled += order/2.0;
  x2scaled += order/2.0; 
  x3scaled += order/2.0;
  x4scaled += order/2.0;
  x5scaled += order/2.0;

  double value = 0.0;
  double d1 = 0.0;
  double d2 = 0.0;
  double d3 = 0.0;
  double d4 = 0.0;
  double d5 = 0.0;

  #pragma omp parallel for collapse(5) reduction(+:value,d1,d2,d3,d4,d5)  
  for(int i = x1grid0; i <= x1grid1; i++){
    for(int j = x2grid0; j <= x2grid1; j++){
      for(int k = x3grid0; k <= x3grid1; k++){
        for(int m = x4grid0; m <= x4grid1; m++){
          for(int n = x5grid0; n <= x5grid1; n++){
            double datapoint = *(data+i*num_x2*num_x3*num_x4*num_x5
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


            value += datapoint*Mn1*Mn2*Mn3*Mn4*Mn5;

            if(f_d){
              double dMn1 = Mn1a - Mn1b;
              double dMn2 = Mn2a - Mn2b;
              double dMn3 = Mn3a - Mn3b;
              double dMn4 = Mn4a - Mn4b;
              double dMn5 = Mn5a - Mn5b;

              d1 += datapoint*dMn1*Mn2*Mn3*Mn4*Mn5
                *(numbin_x1/(max_x1-min_x1));

              d2 += datapoint*Mn1*dMn2*Mn3*Mn4*Mn5
                *(numbin_x2/(max_x2-min_x2));

              d3 += datapoint*Mn1*Mn2*dMn3*Mn4*Mn5
                *(numbin_x3/(max_x3-min_x3));

              d4 += datapoint*Mn1*Mn2*Mn3*dMn4*Mn5 
                *(numbin_x4/(max_x4-min_x4));

              d5 += datapoint*Mn1*Mn2*Mn3*Mn4*dMn5 
                *(numbin_x5/(max_x5-min_x5));
            }
          }
        }
      }
    }
  }
  
  r_value = value;
  r_d1 = d1;
  r_d2 = d2;
  r_d3 = d3;
  r_d4 = d4;
  r_d5 = d5;
}
