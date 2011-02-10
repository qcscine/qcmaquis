#ifndef AMBIENT_CORE_COHERENCY_H
#define AMBIENT_CORE_COHERENCY_H
#include <list>
#include <map>
#include <vector>


namespace ambient { namespace core {

    class coherency_table_entry {


    };

    class coherency_table {
    private:
        coherency_table();
        coherency_table(coherency_table const&);
        coherency_table& operator=(coherency_table const&);
    public:
        ~coherency_table();
        static coherency_table& instance();
        std::map< unsigned long long int, std::vector<coherency_table_entry> > contents;
//      unsigned long long int can be updated to smth in future if not enought
        void update(unsigned long long int id, std::vector<coherency_table_entry> content);

    };

} }
#endif
