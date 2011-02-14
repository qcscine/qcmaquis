#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"

namespace ambient{ namespace core{

    layout_table_entry::layout_table_entry(int owner, int i, int j, int k)
    : owner(owner), i(i), j(j), k(k) { }

    layout_table::~layout_table(){} // don't forget to delete table entries here
    layout_table::layout_table(void_pt_s* object) 
    : object(object)
    {
        if(ambient::is_master()){
            this->reserved_x = 0;
            this->reserved_y = 0;
            update_map();
        }
    }

    void layout_table::update_map(std::vector<layout_table_entry>* update){
        if(update == NULL){
            int y_size = this->object->dim.y / (this->object->group_dim().y*this->object->item_dim().y);
            int x_size = this->object->dim.x / (this->object->group_dim().x*this->object->item_dim().x);
            if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
            for(int i = 0; i < y_size; i++){
                if(i >= this->reserved_y) map.push_back(std::vector<layout_table_entry>());
                for(int j = 0; j < x_size; j++){
                    if(j >= this->reserved_x || i >= this->reserved_y) 
                        map[i].push_back(layout_table_entry(-1, i, j));
                }
            }
            if(x_size > this->reserved_x) this->reserved_x = x_size;
            if(y_size > this->reserved_y) this->reserved_y = y_size;
        }else{
// run through update

        }
    }

    layout_table_entry& layout_table::operator()(const int i, const int j, const int k){
// donno if needed...
        if(ambient::is_master())
            return map[i][j];
        else
            printf("() operation not defined for slave procs\n");
    }

    void layout_table::add_segment_entry(int owner, const int i, const int j, const int k){
        segment.push_back(layout_table_entry(owner, i, j, k));
    }


} }
