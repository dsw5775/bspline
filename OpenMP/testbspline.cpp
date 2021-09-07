#include <iostream>
#include "bspline.h"
#include "common.h"

using std::cout;

int main(int argc, char *argv[])
/*
  testing bspline_1d to bspline_5d with order 6
  you can change order with command line parameter, e.g. ./bspline 4
*/
{
  double start, end;
  
  int order = 6;
  if(argc > 1 && isdigit(argv[1][0])) order = atol(argv[1]);

  int num_x1 = 100+1;
  double min_x1 = 0.0, max_x1 = 1.0;
  double width_x1 = (max_x1-min_x1) / (num_x1-1);
  double x1 = 0.5;

  double value, d1, d2;

  // test bspline_1d for data generated from function x^2
  double *data1 = new double[num_x1];
  for(int i = 0; i < num_x1; i++)
    data1[i] = (width_x1*i) * (width_x1*i);

  start = seconds();
  bspline_1d(data1, min_x1, max_x1,
    num_x1, x1, order, value, 
    d1, true, d2, true);
  end = seconds();
  
  cout << std::fixed;
  cout << "bspline_1d test:\n";
  cout << "Calculated:\n";
  cout << value << '\t' << d1 << '\t' << d2 << '\n';
  cout << "Expected:\n";
  cout << x1*x1 << '\t' << 2*x1 << '\t' << 2.0 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";
  delete[] data1;

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
  double *data2 = new double[num_x1*num_x2];
  double x2 = 4000.0;

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      data2[i*num_x2+j] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
        (min_x2+width_x2*j);

  start = seconds();
  bspline_2d(data2, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    order, value, d1, d2, true);
  end = seconds();
  
  cout << std::fixed;
  cout << "bspline_2d test:\n";
  cout << "Calculated:\n";
  cout << value << '\t' << d1 << '\t' << d2 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2 << '\t' << 2*x1*x2 << '\t' << x1*x1 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";
  delete[] data2;

  /*
    test bspline_3d for data generated from function x1^2*x2*x3
  */
  int num_x3 = 24+1;
  double min_x3 = -0.01, max_x3 = 0.01;
  double width_x3 = (max_x3-min_x3) / (num_x3-1);
  double x3 = 0.005;
  double d3 = 0.0;
  double *data3 = new double[num_x1*num_x2*num_x3];

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        data3[i*num_x2*num_x3+j*num_x3+k] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
          (min_x2+width_x2*j) * (min_x3+width_x3*k);

  start = seconds();
  bspline_3d(data3, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    order, value, d1, d2, d3, true);
  end = seconds();
  
  cout << std::fixed;
  cout << "bspline_3d test:\n";
  cout << "Calculated:\n";
  cout << value << '\t' << d1 << '\t' << d2 << '\t' << d3 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3 << '\t' << 2*x1*x2*x3 << '\t' << x1*x1*x3 << '\t' << x1*x1*x2 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';      
  cout << "\n\n";
  delete[] data3;

  /*
    test bspline_4d for data generated from function x1^2*x2*x3*x4
  */
  int num_x4 = 24+1;
  double min_x4 = -10000.0, max_x4 = 10000.0;
  double width_x4 = (max_x4-min_x4) / (num_x4-1);
  double x4 = 6000.0;
  double d4 = 0.0;
  double *data4 = new double[num_x1*num_x2*num_x3*num_x4];

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        for(int m = 0; m < num_x4; m++)
          data4[i*num_x2*num_x3*num_x4+j*num_x3*num_x4+k*num_x4+m] = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
            (min_x2+width_x2*j) * (min_x3+width_x3*k) * (min_x4+width_x4*m);

  start = seconds();
  bspline_4d(data4, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    order, value, d1, d2, d3, d4, true);
  end = seconds();
  
  cout << std::fixed;
  cout << "bspline_4d test:\n";
  cout << "Calculated:\n";
  cout << value << '\t' << d1 << '\t' << d2 << '\t' << d3 << '\t' << d4 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3*x4 << '\t' << 2*x1*x2*x3*x4 << '\t' << x1*x1*x3*x4 
    << '\t' << x1*x1*x2*x4 << '\t' << x1*x1*x2*x3 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';      
  cout << "\n\n";
  delete[] data4;

  /*
    test bspline_5d for data generated from function x1^2*x2*x3*x4*x5
  */
  int num_x5 = 24+1;
  double min_x5 = -0.01, max_x5 = 0.01;
  double width_x5 = (max_x5-min_x5) / (num_x5-1);
  double x5 = 0.007;
  double d5 = 0.0;
  double *data5 = new double[num_x1*num_x2*num_x3*num_x4*num_x5];

  for(int i = 0; i < num_x1; i++)
    for(int j = 0; j < num_x2; j++)
      for(int k = 0; k < num_x3; k++)
        for(int m = 0; m < num_x4; m++)
          for(int n = 0; n < num_x5; n++)
            data5[i*num_x2*num_x3*num_x4*num_x5+j*num_x3*num_x4*num_x5+k*num_x4*num_x5+m*num_x5+n] 
              = (min_x1+width_x1*i) * (min_x1+width_x1*i) * 
              (min_x2+width_x2*j) * (min_x3+width_x3*k) * (min_x4+width_x4*m) * (min_x5+width_x5*n);

  start = seconds();
  bspline_5d(data5, min_x1, max_x1, num_x1, x1,
    min_x2, max_x2, num_x2, x2,
    min_x3, max_x3, num_x3, x3,
    min_x4, max_x4, num_x4, x4,
    min_x5, max_x5, num_x5, x5,
    order, value, d1, d2, d3, d4, d5, true);
  end = seconds();

  cout << std::fixed;
  cout << "bspline_5d test:\n";
  cout << "Calculated:\n";
  cout << value << '\t' << d1 << '\t' << d2 << '\t' << d3 << '\t' << d4 << '\t' << d5 << '\n';
  cout << "Expected:\n";
  cout << x1*x1*x2*x3*x4*x5 << '\t' << 2*x1*x2*x3*x4*x5 << '\t' << x1*x1*x3*x4*x5 
    << '\t' << x1*x1*x2*x4*x5 << '\t' << x1*x1*x2*x3*x5 << '\t' << x1*x1*x2*x3*x4 << '\n';
  cout << "Time: " << end - start << " sec" << '\n';
  cout << "\n\n";
  delete[] data5;

  return 0;
}
