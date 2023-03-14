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
    valuesSort[i] = values[sortVec[i]];
  }
  
  for (i = 0; i < n; i++){
    values[i] = valuesSort[i];
  }
}
void sortStateVars(State &s, std::vector<int> sortVec);
