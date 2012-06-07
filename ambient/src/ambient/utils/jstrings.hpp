#ifndef AMBIENT_UTILS_JSTRINGS
#define AMBIENT_UTILS_JSTRINGS
#define pinned ambient::models::ambient_pin* ,
#define ctxt_select(...) ctxt_select(std::string(std::string() + __VA_ARGS__).c_str());
#define MAX_NUM_CHAR_LEN 10

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

#endif
