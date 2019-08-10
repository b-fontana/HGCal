#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <map>
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"

#include "../interface/types.h"

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
    T tmp = a.at(size-1);
    typename std::vector<T>::iterator it;
    for(it = a.end()-1; it!=a.begin(); --it) {
      int i = it - a.begin();
      a.at(i) = a.at(i-1);
    }
    a.at(0) = tmp;
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
  
template<typename V>
void check_key(const std::string& k, const mapstr<V>& m) {
  if (m.find(k) == m.end()) {
    std::cout << "Inserted key: " << k << std::endl;
    throw std::invalid_argument("The key does not exist.");
  }
}

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }
  
  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); //ensures that start and end are exactly the same as the input
  return linspaced;
}

TGraphAsymmErrors* build_median_profile(TH2D* h);
std::string etastr(std::string);
#endif //UTILS_H
