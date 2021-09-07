#include <iostream>
//#include <cuda_runtime.h>
#include "bspline.h"
#include "common.h"

using std::cout;

int main(int argc, char **argv)
/*
  testing bspline_1d to bspline_5d with order 6
  you can change order with command line parameter, e.g. ./bspline 4
*/
{
  double *value,*d1,*d2,*d3,*d4,*d5;
  CHECK(cudaMallocManaged((void **)&value, sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d1, sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d2, sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d3, sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d4, sizeof(double)));
  CHECK(cudaMallocManaged((void **)&d5, sizeof(double)));

  double r_value, r_d1, r_d2, r_d3, r_d4, r_d5;

  double start, end;
  
  int order = 6;
  if(argc > 1 && isdigit(argv[1][0])) order = atol(argv[1]);

  int num_x1 = 100+1;
  double min_x1 = 0.0, max_x1 = 1.0;
  double width_x1 = (max_x1-min_x1) / (num_x1-1);
  double x1 = 0.5;

  // test bspline_1d for data generated from function x^2
  double *data1;
  CHECK(cudaMallocManaged((void **)&data1, num_x1*sizeof(double)));

  for(int i = 0; i < num_x1; i++)
    data1[i] = (width_x1*i) * (width_x1*i);

  start = seconds();
  bspline_1d<<<1, 1>>>(data1, min_x1, max_x1,
    num_x1, x1, order, value, 
    d1, true, d2, true);
  cudaDeviceSynchronize();
  end = seconds();

  cout << std::fixed;
  cout << "bspline_1d test:\n";
  cout << "Calculated:\n";
  cout << *value << '\t' << *d1 << '\t' << *d2 << '\n';
  cout << "Expected:\n";
  cout << x1*x1 << '\t' << 2*x1 << '\t' << 2.0 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";

  /*
    test bspline_2d for data generated from function x1^2*x2
  */
  num_x1 = 100; 
  min_x1 = 0.005;
  max_x1 = 0.995;
  width_x1 = (max_x1-min_x1) / (num_x1-1);
  x1 = 0.3;
  int num_x2 = 24+1;
  double min_x2 = -10000.0, max_x2 = 10000.0;
  double width_x2 = (max_x2-min_x2) / (num_x2-1);
  double *data2; // new double[num_x1*num_x2];
  double x2 = 4000.0;
  CHECK(cudaMallocManaged((void **)&data2, num_x1*num_x2*sizeof(double)));

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      data2[i*num_x2+j] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
        (min_x2+width_x2*j);

  start = seconds();
  bspline_2d<<<1, 1>>>(data2, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    order, value, d1, d2, true);
  cudaDeviceSynchronize();
  end = seconds();

  cout << std::fixed;
  cout << "bspline_2d test:\n";
  cout << "Calculated:\n";
  cout << *value << '\t' << *d1 << '\t' << *d2 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2 << '\t' << 2*x1*x2 << '\t' << x1*x1 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";

  /*
    test bspline_3d for data generated from function x1^2*x2*x3
  */
  int num_x3 = 24+1;
  double min_x3 = -0.01, max_x3 = 0.01;
  double width_x3 = (max_x3-min_x3) / (num_x3-1);
  double x3 = 0.005;
  double *data3; // new double[num_x1*num_x2*num_x3];
  CHECK(cudaMallocManaged((void **)&data3, num_x1*num_x2*num_x3*sizeof(double)));

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        data3[i*num_x2*num_x3+j*num_x3+k] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
          (min_x2+width_x2*j) * (min_x3+width_x3*k);

  start = seconds();
  /*
  bspline_3d<<<1, order*order*order>>>(data3, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    order, value, d1, d2, d3, true);
  cudaDeviceSynchronize();
  */
  bspline_3d_ex(data3, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    order, r_value, r_d1, r_d2, r_d3, true);
  end = seconds();

  cout << std::fixed;
  cout << "bspline_3d test:\n";
  cout << "Calculated:\n";
  cout << r_value << '\t' << r_d1 << '\t' << r_d2 << '\t' << r_d3 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3 << '\t' << 2*x1*x2*x3 << '\t' << x1*x1*x3 << '\t' << x1*x1*x2 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";

  /*
    test bspline_4d for data generated from function x1^2*x2*x3*x4
  */
  int num_x4 = 24+1;
  double min_x4 = -10000.0, max_x4 = 10000.0;
  double width_x4 = (max_x4-min_x4) / (num_x4-1);
  double x4 = 6000.0;
  double *data4; // new double[num_x1*num_x2*num_x3*num_x4];
  CHECK(cudaMallocManaged((void **)&data4, num_x1*num_x2*num_x3*num_x4*sizeof(double)));

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        for(int m = 0; m < num_x4; m++)
          data4[i*num_x2*num_x3*num_x4+j*num_x3*num_x4+k*num_x4+m] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
            (min_x2+width_x2*j) * (min_x3+width_x3*k) * (min_x4+width_x4*m);

  start = seconds();
  bspline_4d_ex(data4, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    order, r_value, r_d1, r_d2, r_d3, r_d4, true);
  end = seconds();

  cout << std::fixed;
  cout << "bspline_4d test:\n";
  cout << "Calculated:\n";
  cout << r_value << '\t' << r_d1 << '\t' << r_d2 << '\t' << r_d3 << '\t' << r_d4 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3*x4 << '\t' << 2*x1*x2*x3*x4 << '\t' << x1*x1*x3*x4 
    << '\t' << x1*x1*x2*x4 << '\t' << x1*x1*x2*x3 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";

  /*
    test bspline_5d for data generated from function x1^2*x2*x3*x4*x5
  */
  int num_x5 = 24+1;
  double min_x5 = -0.01, max_x5 = 0.01;
  double width_x5 = (max_x5-min_x5) / (num_x5-1);
  double x5 = 0.007;

  double *data5;
  int nBytes = num_x1*num_x2*num_x3*num_x4*num_x5*sizeof(double);
  CHECK(cudaMallocManaged((void **)&data5, nBytes));
  
  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        for(int m = 0; m < num_x4; m++)
          for(int n = 0; n < num_x5; n++)
            data5[i*num_x2*num_x3*num_x4*num_x5+j*num_x3*num_x4*num_x5+k*num_x4*num_x5+m*num_x5+n] 
              = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
              (min_x2+width_x2*j) * (min_x3+width_x3*k) * (min_x4+width_x4*m) * (min_x5+width_x5*n);

  start = seconds();
  bspline_5d_ex(data5, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    min_x5, max_x5, num_x5, x5,
    order, r_value, r_d1, r_d2, r_d3, r_d4, r_d5, true);
  end = seconds();

  cout << std::fixed;
  cout << "bspline_5d test:\n";
  cout << "Calculated:\n";
  cout << r_value << '\t' << r_d1 << '\t' << r_d2 << '\t' << r_d3 
    << '\t' << r_d4 << '\t' << r_d5 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3*x4*x5 << '\t' << 2*x1*x2*x3*x4*x5 << '\t' << x1*x1*x3*x4*x5 
    << '\t' << x1*x1*x2*x4*x5 << '\t' << x1*x1*x2*x3*x5 << '\t' << x1*x1*x2*x3*x4 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";

  // free device global memory
  CHECK(cudaFree(data1));
  CHECK(cudaFree(data2));
  CHECK(cudaFree(data3));
  CHECK(cudaFree(data4));
  CHECK(cudaFree(data5));
  CHECK(cudaFree(value));
  CHECK(cudaFree(d1));
  CHECK(cudaFree(d2));
  CHECK(cudaFree(d3));
  CHECK(cudaFree(d4));
  CHECK(cudaFree(d5));

  // reset device
  CHECK(cudaDeviceReset());

  return 0;
}
