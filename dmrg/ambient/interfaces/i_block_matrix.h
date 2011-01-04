#ifndef AMBIENT_BLOCK_MATRIX_IFACE
#define AMBIENT_BLOCK_MATRIX_IFACE

class i_block_matrix{
   public:
   virtual void* memory_pointer() = 0;
//   virtual void* block_memory_pointer(int k) = 0;
//   virtual size_t block_lda(int k) = 0;
//   virtual bool is_coherent() = 0;
//   virtual void set_coherency(bool coherent) = 0;
//   virtual int num_blocks() = 0;
//   virtual size_t block_size(int k) = 0;
//   virtual breakdown* breakdown_set() = 0;
   virtual ~i_block_matrix(){}
};

#endif
