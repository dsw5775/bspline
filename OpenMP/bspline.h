#pragma once
#include <stdio.h>
#include <stdlib.h>

void bspline_1d(double *data, double min_x1, double max_x1,
  int num_x1, double x1, int order, double& r_value, 
  double& r_d1, bool f_d1, double& r_d2, bool f_d2);
  
void bspline_2d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  int order, double& r_value, double& r_d1, double& r_d2, 
  bool f_d);
    
void bspline_3d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3,
  bool f_d);

void bspline_4d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4,
  bool f_d);  

void bspline_5d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4, double& r_d5,
  bool f_d);  

inline double Mn(int n, double u)
/* nth-order B-spline */
{
  if(n < 2) {
    printf("Order of B-spline is less than 2\n");
    exit(-1);
  }

  if(n == 2) {
    if(u < 0.0 || u > 2.0) return 0.0;

    if(u < 1.0)
      return u;
    else
      return 2.0 - u;
  }

  return u/(n-1)*Mn(n-1, u) + (n-u)/(n-1)*Mn(n-1, u-1);
}

inline double dMn(int n, double u)
/* first derivative */
{
  /* 
  // taken care of by Mn 
  if(n < 3) {
    printf("Order of B-spline is less than 3 for first derivative\n");
    exit(-1);
  }
  */

  return Mn(n-1, u) - Mn(n-1, u-1);
}

inline double dMn2(int n, double u)
/* second derivative */
{
  /* 
  // taken care of by Mn
  if(n < 4) {
    printf("Order of B-spline is less than 4 for second derivative\n");
    exit(-1);
  }
  */

  return Mn(n-2, u) - 2*Mn(n-2, u-1) + Mn(n-2, u-2);
}
