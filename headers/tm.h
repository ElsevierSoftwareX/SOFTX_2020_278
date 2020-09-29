#ifndef __TM_H
#define __TM_H

#include <sstream>
#include <vector>

using t1D = std::vector<double>;
using t2D = std::vector<t1D>;
using t3D = std::vector<t2D>;

template <class T> T sqr(T x)
{
  return x*x;
}
//-----------------------------------
template <typename T>
std::string Str(const T & t) {
   std::ostringstream os;
   os << t;
   return os.str();
}

#endif // TM_H
