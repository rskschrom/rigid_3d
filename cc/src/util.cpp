#include <iostream>
#include <math.h>
#include <vector>
#include "state.h"

// sort array via sort vector
template <typename T>
void sortWithVector(std::vector<T> &values, std::vector<int> sortVec){
  int i;
  int n = sortVec.size();
  std::vector<T> valuesSort(n);
  
  for (i = 0; i < n; i++){
    valuesSort[sortVec[i]] = values[i];    
  }
  
  for (i = 0; i < n; i++){
    values[i] = valuesSort[i];
  }
}

// sort all of the state variable vectors
void sortStateVars(State &s, std::vector<int> sortVec){
  std::vector<float> x = s.x;
  std::vector<float> y = s.y;
  std::vector<int> hi = s.hi;
  
  sortWithVector<float>(x, sortVec);
  sortWithVector<float>(y, sortVec);
  sortWithVector<int>(hi, sortVec);
  
  s.x = x;
  s.y = y;
  s.hi = hi;
}
