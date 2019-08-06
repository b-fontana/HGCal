#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <unordered_map>
#include <cstddef>

using char_ = char;
using bool_ = bool;
using float_ = float;
using int_ = int;
using uint_ = unsigned int;

template <typename T> 
using vec1d = std::vector<T>;

template <typename T> 
using vec2d = std::vector< std::vector<T> >;

template <typename T> 
using vec3d = std::vector< std::vector< std::vector<T> > >;

template <typename T>
class vector3d {
 public:
 vector3d(std::size_t d1=0, std::size_t d2=0, std::size_t d3=0, T const& t=T()):
  d1(d1), d2(d2), d3(d3), data(d1*d2*d3, t) {}
  
  T& operator()(std::size_t i, std::size_t j, std::size_t k) {
    return data[i*d2*d3 + j*d3 + k];
  }
  T const& operator()(std::size_t i, std::size_t j, std::size_t k) const {
    return data[i*d2*d3 + j*d3 + k];
  }
  
 private:
  std::size_t d1, d2, d3;
  vec1d<T> data;
};

template <typename T>
using mapstr = std::unordered_map<std::string, T>;

#endif //TYPES_H
