#include "ambient/ambient.h"
#include "ambient/interface/p_profile.h"

namespace ambient {

    p_profile_s::p_profile_s():reserved_x(0),reserved_y(0),group_id(0),id(0){};

    p_profile* p_profile_s::dereference(){
        while((this->profile = this->profile->profile) != this->profile->profile);
        return this->profile; // todo - deallocate proxy objects
    }

    void p_profile_s::set_id(std::pair<unsigned int*,size_t> group_id){
        this->group_id = group_id.first;
        this->layout = new core::layout_table();
        this->id = void_pt_map.insert(group_id.first, group_id.second, this->layout);
    }

    p_profile & p_profile_s::operator>>(dim3 dim_distr) 
    {
        this->specific = true;
        this->dim_distr = dim_distr;
        this->dim_group = NULL;
        this->dim_gpu = NULL;
        return *(p_profile*)this;
    }
    p_profile & p_profile_s::operator,(dim3 dim) 
    {
        if(this->dim_group == NULL){
            this->dim_group = dim;
            this->regroup();
        }else if(this->dim_gpu == NULL){
            this->dim_gpu = dim;
        }
        return *(p_profile*)this;
    }

    p_profile & operator>>(p_profile* instance, dim3 dim_distr) {
        return *instance >> dim_distr;
    }

    void p_profile_s::regroup(){
        if(!this->proxy){
            int y_size = this->dim.y / (this->group_dim().y*this->item_dim().y);
            int x_size = this->dim.x / (this->group_dim().x*this->item_dim().x);
            if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
            for(int i = 0; i < y_size; i++){
                if(i >= this->reserved_y) skeleton.push_back(std::vector<workgroup*>());
                for(int j = 0; j < x_size; j++){
                    if(j >= this->reserved_x || i >= this->reserved_y) 
                        skeleton[i].push_back(new workgroup(&profile, i, j));
                }
            }
            if(x_size > this->reserved_x) this->reserved_x = x_size;
            if(y_size > this->reserved_y) this->reserved_y = y_size;
        }
    }

    workgroup* p_profile_s::group(int i, int j, int k){
        if(this->proxy){
            return new workgroup(&profile, i, j, k);
        }else{
            int x_size = this->dim.x / (this->group_dim().x*this->item_dim().x);
            int y_size = this->dim.y / (this->group_dim().y*this->item_dim().y);
            int z_size = this->dim.z / (this->group_dim().z*this->item_dim().z);

            if(i >= y_size || j >= x_size || k >= z_size) printf("Warning: accessing group that is out of range\n");
            return this->skeleton[i][j];
        }
    }

    dim3 p_profile_s::grid_dim(){
        int x_size = this->dim.x / (this->group_dim().x * this->item_dim().x);
        int y_size = this->dim.y / (this->group_dim().y * this->item_dim().y);
        int z_size = this->dim.z / (this->group_dim().z * this->item_dim().z);
        return dim3(x_size, y_size, z_size);
    }

    dim3 p_profile_s::group_dim(){
        if(this->specific) return this->dim_group;
        else return engine.group_dim();
    }

    dim3 p_profile_s::item_dim(){
        return engine.item_dim();
    }

}
