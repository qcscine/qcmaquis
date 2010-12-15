#include "utils/interface_macros.h"

class i_block_matrix{
   interface:
   virtual void* memory_pointer() = 0;

   virtual ~i_block_matrix(){}
};

