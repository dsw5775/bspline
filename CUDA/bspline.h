#pragma once

#include <assert.h>

#define FULL_MASK 0xffffffff
#define BLOCK_SIZE 512 // be a multiple of warpSize for not to cause inactive threads in shuffle 

void bspline_3d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3,
  bool f_d);

void bspline_4d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4,
  bool f_d);

void bspline_5d_ex(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double& r_value, double& r_d1, double& r_d2, double& r_d3, double& r_d4, double& r_d5,
  bool f_d);

__global__ void bspline_1d(double *data, double min_x1, double max_x1,
  int num_x1, double x1, int order, double* o_value, 
  double* o_d1, bool f_d1, double* o_d2, bool f_d2);
  
__global__ void bspline_2d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  int order, double* o_value, double* o_d1, double* o_d2,
  bool f_d);
    
__global__ void bspline_3d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  int order, double* o_value, double* o_d1, double* o_d2, double* o_d3,
  bool f_d);

__global__ void bspline_4d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  int order, double* value, double* d1, double* d2, double* d3, double* d4,
  bool f_d);

__global__ void bspline_5d(double *data,
  double min_x1, double max_x1, int num_x1, double x1, 
  double min_x2, double max_x2, int num_x2, double x2, 
  double min_x3, double max_x3, int num_x3, double x3,
  double min_x4, double max_x4, int num_x4, double x4,
  double min_x5, double max_x5, int num_x5, double x5,  
  int order, double* o_value, double* o_d1, double* o_d2, double* o_d3, double* o_d4, double* o_d5,
  bool f_d);

__inline__ __device__
double Mn(int n, double u)
/* nth-order B-spline */
{
  /*
    without __trap() or assert, inline (or even noinline) recursive won't work
  */

  if(n < 2){
    printf("Order of B-spline is less than 2\n");
    __trap();
  } 
  
  //assert(n >= 2);

  if(n == 2) {
    if(u < 0.0 || u > 2.0) return 0.0;

    if(u < 1.0)
      return u;
    else
      return 2.0 - u;
  }

  return u/(n-1)*Mn(n-1, u) + (n-u)/(n-1)*Mn(n-1, u-1);
}

__inline__ __device__
double dMn(int n, double u)
/* first derivative */
{
  return Mn(n-1, u) - Mn(n-1, u-1);
}

__inline__ __device__
double dMn2(int n, double u)
/* second derivative */
{
  return Mn(n-2, u) - 2*Mn(n-2, u-1) + Mn(n-2, u-2);
}

__inline__ __device__ double warpReduceSum(double localSum)
{
  localSum += __shfl_down_sync(FULL_MASK, localSum, 16);
  localSum += __shfl_down_sync(FULL_MASK, localSum, 8);
  localSum += __shfl_down_sync(FULL_MASK, localSum, 4);
  localSum += __shfl_down_sync(FULL_MASK, localSum, 2);
  localSum += __shfl_down_sync(FULL_MASK, localSum, 1);

  return localSum;
}

__inline__ __device__
double blockReduceSum(double val) {

  static __shared__ double shared[32]; // Shared mem for 32 partial sums
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  val = warpReduceSum(val);     // Each warp performs partial reduction

  if (lane==0) shared[wid]=val; // Write reduced value to shared memory

  __syncthreads();              // Wait for all partial reductions

  //read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;

  if (wid==0) val = warpReduceSum(val); //Final reduce within first warp

  return val;
}

