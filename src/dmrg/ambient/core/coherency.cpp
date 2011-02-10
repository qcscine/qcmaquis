#include "ambient/core/coherency.h"

namespace ambient{ namespace core{

    coherency_table& coherency_table::instance()
    {
        static coherency_table* singleton = NULL;
        if(!singleton) singleton = new coherency_table();
        return *singleton;
    }

    coherency_table::coherency_table(){}

    coherency_table::~coherency_table(){} // don't forget to delete table entries here

    void update(unsigned long long int id, std::vector<coherency_table_entry> content){

// now let's find this id in our map table

// now let's update it

    }


} }
