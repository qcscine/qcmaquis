#ifndef AMBIENT_CORE_LAYOUT_H
#define AMBIENT_CORE_LAYOUT_H
#include <list>
#include <map>
#include <vector>

namespace ambient { 
    class p_profile;
    class p_profile_s;

namespace core {
    typedef ambient::p_profile   void_pt;
    typedef ambient::p_profile_s void_pt_s;

    class layout_table_entry {
    public:
        layout_table_entry(int owner, int i, int j = 0, int k = 0);
        int i;
        int j;
        int k;
        int owner;
    };

    class layout_table {
    public:
        layout_table(void_pt_s* object);
        ~layout_table();

        void update_map(std::vector<layout_table_entry>* update = NULL);
        layout_table_entry* operator()(const int i, const int j = 0, const int k = 0);

        void add_segment_entry(int owner, int i, int j = 0, int k = 0);
        void update_map_entry(int owner, int i, int j = 0, int k = 0);

        void_pt_s* object;
        std::vector< std::vector<layout_table_entry*> > map;
        size_t reserved_x;
        size_t reserved_y;
        std::vector<layout_table_entry> segment;
        size_t count;
    };

} }
#endif
