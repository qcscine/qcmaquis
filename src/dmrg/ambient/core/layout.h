#ifndef AMBIENT_CORE_LAYOUT_H
#define AMBIENT_CORE_LAYOUT_H
#include <list>
#include <map>
#include <vector>


namespace ambient { namespace core {

    class layout_table_entry {


    };

    class layout_table {
    public:
        layout_table();
        ~layout_table();
        std::vector<layout_table_entry> contents;
//      unsigned long long int can be updated to smth in future if not enought
        void update(std::vector<layout_table_entry> content);
    };

} }
#endif
