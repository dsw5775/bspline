#pragma once
#include <stdio.h>
#include <stdlib.h>

template <typename T>
void bspline_1d(T *data, T min_x1, T max_x1,
  int num_x1, T x1, int order, T& r_value, 
  T& r_d1, bool f_d1, T& r_d2, bool f_d2);
  
template <typename T>
void bspline_2d(T *data,
  T min_x1, T max_x1, int num_x1, T x1, 
  T min_x2, T max_x2, int num_x2, T x2, 
  int order, T& r_value, T& r_d1, T& r_d2, 
  bool f_d);
    
template <typename T>
void bspline_3d(T *data,
  T min_x1, T max_x1, int num_x1, T x1, 
  T min_x2, T max_x2, int num_x2, T x2, 
  T min_x3, T max_x3, int num_x3, T x3,
  int order, T& r_value, T& r_d1, T& r_d2, T& r_d3,
  bool f_d);

template <typename T>
void bspline_4d(T *data,
  T min_x1, T max_x1, int num_x1, T x1, 
  T min_x2, T max_x2, int num_x2, T x2, 
  T min_x3, T max_x3, int num_x3, T x3,
  T min_x4, T max_x4, int num_x4, T x4,
  int order, T& r_value, T& r_d1, T& r_d2, T& r_d3, T& r_d4,
  bool f_d);  

template <typename T>
void bspline_5d(T *data,
  T min_x1, T max_x1, int num_x1, T x1, 
  T min_x2, T max_x2, int num_x2, T x2, 
  T min_x3, T max_x3, int num_x3, T x3,
  T min_x4, T max_x4, int num_x4, T x4,
  T min_x5, T max_x5, int num_x5, T x5,  
  int order, T& r_value, T& r_d1, T& r_d2, T& r_d3, T& r_d4, T& r_d5,
  bool f_d);  

template <typename T>
inline T Mn(int n, T u)
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

template <typename T>
inline T dMn(int n, T u)
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

template <typename T>
inline T dMn2(int n, T u)
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

#include "bspline.ipp"
