#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <tuple>

#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"

/*shifts vector one value to the right*/
template <typename T> 
std::vector<T> shift(std::vector<T> a) {
  unsigned int size = a.size();
  T tmp = a[size-1];
  typename std::vector<T>::iterator it;
  for(it = a.end(); it!=a.begin(); --it) {
    int i = it - a.begin(); 
    a[i] = a[i-1];
  }
  a[0] = tmp;
  return a;
}

TGraphAsymmErrors* build_median_profile(TH2D* h);

#endif //UTILS_H
