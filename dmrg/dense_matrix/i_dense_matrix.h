#include "utils/interface_macros.h"

class i_dense_matrix{
   interface:
   virtual void* memory_pointer() = 0;

   virtual ~i_dense_matrix(){}
};
