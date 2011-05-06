/**
  Some common macros
*/

#define NUM 32;
//#define NUM 1;

// BASE 255

#define LOG_BASE			 0x8   // 8 
#define LOG_BASE_HALF		 0x4   // 4
#define BASE				 0x100 // 256
#define BASE_HALF			 0x10  // 16
#define BASE_MINUS2			 0xFE  // 254   
#define MINUS_BASE_PLUS2	 -254
#define MASK_DOWN			 0xF
#define MASK_UP				 0xF0

/*
//BASE 2**30
#define LOG_BASE			 30   
#define LOG_BASE_HALF		 15   
#define BASE				 0x7FFFFFFF 
#define BASE_HALF			 0x3FFFFFFF 
#define BASE_MINUS2			 0x7FFFFFFD
#define MINUS_BASE_PLUS2	 -0x7FFFFFFD
#define MASK_DOWN			 0xFFFF
#define MASK_UP				 0xFFFF0000
*/

namespace vli {

    typedef std::size_t size_int;
    typedef std::size_t size_type;
}
