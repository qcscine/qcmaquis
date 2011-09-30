#ifndef MAQUIS_ZOUT_H
#define MAQUIS_ZOUT_H


#ifdef _AMBIENT
  #include <ambient/utils/zout.hpp>
  #define zout ambient::zout
#else
  #include <iostream>
  #define zout std::cout
#endif


#endif
