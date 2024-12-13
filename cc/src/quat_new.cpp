#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

// mutliplication of quaternions
void multiply(std::vector<float> q1, std::vector<float> q2, std::vector<float> &q3){
  float r1 = q1[0], v1x = q1[1], v1y = q1[2], v1z = q1[3];
  float r2 = q2[0], v2x = q2[1], v2y = q2[2], v2z = q2[3];
  
  q3[0] = r1*r2-v1x*v2x-v1y*v2y-v1z*v2z;
  q3[1] = r1*v2x+r2*v1x+v1y*v2z-v1z*v2y;
  q3[2] = r1*v2y+r2*v1y-v1x*v2z+v1z*v2x;
  q3[3] = r1*v2z+r2*v1z+v1x*v2y-v1y*v2x;
}

// conjugate of quaternion
void conj(std::vector<float> q, std::vector<float> &cq){
  cq[0] = q[0];
  cq[1] = -q[1];
  cq[2] = -q[2];
  cq[3] = -q[3];
}

// rotate pure quaternion (vector) with unit quaternion (versor)
void vect_rotate(std::vector<float> v, std::vector<float> q, std::vector &vr){
  float a, b, c;

  a = q[0]*q[0]-
      q[1]*q[1]-
      q[2]*q[2]-
      q[3]*q[3];  
  b = 2.*(v[0]*q[1]+
          v[1]*q[2]+
          v[2]*q[3]);
  c = -2.*q[0];
  
  // vr = a*p+b*v+c*pxv
  vr[0] = a*v[0]+b*q[1]+c*(v[1]*q[3]-v[2]*q[2]);
  vr[1] = a*v[1]+b*q[2]+c*(v[2]*q[1]-v[0]*q[3]);
  vr[2] = a*v[2]+b*q[3]+c*(v[0]*q[2]-v[1]*q[1]);
}
