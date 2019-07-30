#ifndef TYPES_H
#define TYPES_H

#include <vector>

#include "TF1.h"

template <typename T> 
using vec1d = std::vector<T>;

template <typename T> 
using vec2d = std::vector< std::vector<T> >;

template <typename T> 
using vec3d = std::vector< std::vector< std::vector<T> > >;

using float_ = float;
using int_ = int;
using uint_ = unsigned int;
using mapfunc = std::map<std::string, TF1* >;

/*
typedef std::vector< std::vector<float> > vec2d;
typedef std::vector< std::vector< std::vector<float> > > vec3d;

typedef std::vector<TGraph*> vec_graph;
typedef std::map< std::string, TF1* > map_func; //calibration storage map 
*/

#endif //TYPES_H
