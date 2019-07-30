#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <map>

#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"

template <typename T>
class VecOps {
 private:
  std::vector<T> a;
  unsigned int s;

 public:
 VecOps(std::vector<T> a_): a(a_) {
    this->s = a.size();
  };
  
  /*shifts vector one value to the right*/  
  std::vector<T> shift() {
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

  std::vector<T> average_array(T last) {
    typename std::vector<T>::iterator it;
    for(it=a.begin(); a.end()-it>1; ++it) {
      int idx = std::distance(a.begin(), it);
      a[idx] = (a[idx]+a[idx+1])/2;
    }
    a[s-1] = last;
    return a;
  }
}; 
  
template<typename V, typename K>
  void check_key(const std::map<K,V>& m, const K& k) {
  if (m.find(k) == m.end()) {
    std::cout << "Inserted key: " << k << std::endl;
    throw std::invalid_argument("The key does not exist.");
  }
}

TGraphAsymmErrors* build_median_profile(TH2D* h);
std::string etastr(std::string);
#endif //UTILS_H
