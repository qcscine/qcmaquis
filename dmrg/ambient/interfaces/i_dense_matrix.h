class i_dense_matrix{
   public:
   virtual void* memory_pointer() = 0;
//   virtual void* block_memory_pointer(int i, int j) = 0;
//   virtual size_t lda() = 0;
//   virtual bool is_coherent() = 0;
//   virtual void set_coherency(bool coherent) = 0;
//   virtual size_t size() = 0;
//   virtual size_t block_size() = 0;
//   virtual breakdown* breakdown_set() = 0;
   virtual ~i_dense_matrix(){}
};

