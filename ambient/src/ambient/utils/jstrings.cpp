#include "ambient/utils/jstrings.h"

namespace ambient { 

    std::string operator+(std::string lhs, double rhs){
        char rhs_str[MAX_NUM_CHAR_LEN];
        if(rhs - (int)rhs) sprintf(rhs_str,"%.1f",rhs);
        else sprintf(rhs_str,"%d",(int)rhs);
        lhs += rhs_str;
        return lhs;
    }

    std::string operator+(std::string lhs, std::pair<size_t*,size_t> rhs){
        return (lhs + (int)(*rhs.first) + "-" + (int)rhs.second);
    }

}
