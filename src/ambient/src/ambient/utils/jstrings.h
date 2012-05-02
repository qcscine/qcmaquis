#ifndef AMBIENT_UTILS_JSTRINGS
#define AMBIENT_UTILS_JSTRINGS
#define pinned ambient::models::ambient_pin* ,
#define ctxt_select(...) ctxt_select(std::string(std::string() + __VA_ARGS__).c_str());
#define MAX_NUM_CHAR_LEN 10

#include <string>

namespace ambient { 
    std::string operator+(std::string lhs, double rhs);
    std::string operator+(std::string lhs, std::pair<size_t*,size_t> rhs);
}

#endif
