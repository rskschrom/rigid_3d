#include <iostream>
#include <math.h>
#include <vector>

std::vector<float> multiply(std::vector<float> q1, std::vector<float> q2);
std::vector<float> multiplyVecQuat(std::vector<float> v, std::vector<float> q);
std::vector<float> conj(std::vector<float> q);
std::vector<float> vecRotate(std::vector<float> v, std::vector<float> q);