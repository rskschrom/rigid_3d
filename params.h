#ifndef Params_H
#define Params_H

//structure for simulation parameters
struct Params{
  int nt = 500;
  float pi = 3.14159265;
  float ds = 0.1;
  float dt = 1.e-2;
  float g = -9.81;
  float fc = 0.5;
};
#endif
