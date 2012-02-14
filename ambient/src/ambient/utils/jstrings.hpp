#define pinned ambient::models::ambient_pin* ,
#define ctxt_select(...) ctxt_select(std::string(std::string() + __VA_ARGS__).c_str());
#define MAX_NUM_CHAR_LEN 10

namespace ambient { 

    // string mangling for queries //
    std::string & operator+(std::string & lhs, double rhs){
        char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
        if(rhs - (int)rhs) sprintf(rhs_str,"%.1f",rhs);
        else sprintf(rhs_str,"%d",(int)rhs);
        lhs += rhs_str;
        free(rhs_str);
        return lhs;
    }

    std::string & operator+(std::string & lhs, int rhs){
        char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
        sprintf(rhs_str,"%d", rhs);
        lhs += rhs_str;
        free(rhs_str);
        return lhs;
    }

    std::string & operator+(const std::string & lhs, double rhs){
        return const_cast<std::string&>(lhs)+rhs;
    }

    std::string & operator+(const std::string & lhs, int rhs){
        return const_cast<std::string&>(lhs)+rhs;
    }

    std::string & operator+(std::string & lhs, const char* rhs){
        return lhs += rhs;
    }

    std::string & operator+(const std::string & lhs, const char* rhs){
        return const_cast<std::string&>(lhs) += rhs;
    }

    std::string & operator+(std::string & lhs, std::pair<size_t*,size_t> rhs){
        lhs + (int)(*rhs.first);
        lhs += "-";
        lhs + (int)(rhs.second);
        return lhs;
    }

    std::string & operator+(const std::string & lhs, std::pair<size_t*,size_t> rhs){
        const_cast<std::string&>(lhs) + (int)(*rhs.first);
        const_cast<std::string&>(lhs) += "-";
        return const_cast<std::string&>(lhs) + (int)rhs.second;
    }

} // namespace ambient
