#ifndef UTILS_H
#define UTILS_H

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

#endif //UTILS_H
