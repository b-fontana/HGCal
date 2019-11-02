#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <exception>

using char_ = char;
using bool_ = bool;
using double_ = double;
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
using vec4d = std::vector< std::vector< std::vector< std::vector<T> > > >;

template <typename T> 
using tup3 = std::tuple<T, T, T>;

template <typename T>
class vector3d {
 public:
 vector3d(std::size_t d1=0, std::size_t d2=0, std::size_t d3=0, T const& t=T()):
  d1(d1), d2(d2), d3(d3), data(d1*d2*d3, t) {}
  
  T& operator()(std::size_t i, std::size_t j, std::size_t k) 
    {
      return data[i*d2*d3 + j*d3 + k];
    }
  T const& operator()(std::size_t i, std::size_t j, std::size_t k) const 
  {
    return data[i*d2*d3 + j*d3 + k];
  }
  T& at(std::size_t i, std::size_t j, std::size_t k)
    {
      if(i>=d1 or j>=d2 or k>=d3)
	throw std::out_of_range("The ("+std::to_string(i)+" "+std::to_string(j)+" "+std::to_string(k)
				+") entry exceeds the allocated vector size.");
      else if(i<0 or j<0 or k<0)
	throw std::out_of_range("The ("+std::to_string(i)+" "+std::to_string(j)+" "+std::to_string(k)
				+") entry must be non-negative.");
      return data[i*d2*d3 + j*d3 + k];
    }

 private:
  std::size_t d1, d2, d3;
  vec1d<T> data;
};

template <typename T> 
using opt = std::optional<T>;

template <typename T, typename U> 
  using umap = std::unordered_map<T,U>;

template <typename T>
using mapstr = umap<std::string, T>;

enum class ParticleType {Photon, Pion, NTypes};
extern umap<ParticleType, std::string> mParticleType;
enum class DetectorRegion {Inner, Central, Outer, NRegions};
extern umap<DetectorRegion, std::string> mDetectorRegion;
enum class Method {ShowerLeakage, BruteForce};
extern umap<Method, std::string> mMethod;

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}

#endif //TYPES_H
