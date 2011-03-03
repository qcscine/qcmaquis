#ifndef AMBIENT_CORE_LAYOUT_H
#define AMBIENT_CORE_LAYOUT_H
#include <list>
#include <map>
#include <vector>

namespace ambient { 
    class p_profile;

namespace core {

    class layout_table_entry {
    public:
        layout_table_entry(); // default constructor
        layout_table_entry(int owner, int i, int j = 0, int k = 0);
        int i;
        int j;
        int k;
        int owner;
        int xowner;
        int get_xowner();
        int get_owner();
    };

    class layout_table {
    public:
        layout_table(p_profile* object);
        ~layout_table();

        void update_map(std::vector<layout_table_entry>* update = NULL);
        layout_table_entry* get_entry(int i, int j = 0, int k = 0);
        layout_table_entry* operator()(int i, int j = 0, int k = 0);

        void add_segment_entry(int owner, int i, int j = 0, int k = 0);
        void update_map_entry(int owner, int i, int j = 0, int k = 0);
        void record(int owner, int i, int j = 0, int k = 0); // general call invoking one above
        void clean();
        void print();
        void apply();

        p_profile* object;
        std::vector< std::vector<layout_table_entry*> > map;
        size_t reserved_x;
        size_t reserved_y;
        std::vector<layout_table_entry> segment;
        size_t count;
        size_t segment_count;
        size_t xsegment_count;
    };

    void apply_change_set(p_profile** profiles, size_t count);
    void perform_forwarding(p_profile** profiles, size_t count);

} }
#endif
