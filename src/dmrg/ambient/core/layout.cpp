#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/auxiliary.hpp"

namespace ambient{ namespace core{

    layout_table_entry::layout_table_entry()
    { }

    layout_table_entry::layout_table_entry(int owner, int i, int j, int k)
    : owner(owner), i(i), j(j), k(k) { }

    layout_table::~layout_table(){} // don't forget to delete table entries here
    layout_table::layout_table(p_profile* object) 
    : object(object), count(0), segment_count(0), xsegment_count(0)
    {
        this->reserved_x = 0;
        this->reserved_y = 0;
        update_map();
    }

    void layout_table::update_map(std::vector<layout_table_entry>* update){
        if(update == NULL){
            int y_size = this->object->dim.y / (this->object->get_group_dim().y*this->object->get_item_dim().y);
            int x_size = this->object->dim.x / (this->object->get_group_dim().x*this->object->get_item_dim().x);
            if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
            for(int i = 0; i < y_size; i++){
                if(i >= this->reserved_y) map.push_back(std::vector<layout_table_entry*>());
                for(int j = 0; j < x_size; j++){
                    if(j >= this->reserved_x || i >= this->reserved_y) 
                        map[i].push_back(NULL); //new layout_table_entry(-1, i, j);
                }
            }
            if(x_size > this->reserved_x) this->reserved_x = x_size;
            if(y_size > this->reserved_y) this->reserved_y = y_size;
        }else{
// run through update

        }
    }

    layout_table_entry* layout_table::get_entry(int i, int j, int k){
        return map[i][j];
    }

    layout_table_entry* layout_table::operator()(int i, int j, int k){
        return map[i][j];
    }

    void layout_table::add_segment_entry(int owner, int i, int j, int k){
        if(segment_count == segment.size()) segment.resize(segment_count+1);
        segment[segment_count++] = layout_table_entry(owner, i, j, k);
    }

    void layout_table::update_map_entry(int owner, int i, int j, int k){
        if(map[i][j] == NULL) this->count++;
        map[i][j] = new layout_table_entry(owner, i, j, k);
    }
    void layout_table::record(int owner, int i, int j, int k) {
        if(ambient::scope.master()){ 
            update_map_entry(owner, i, j, k);
            add_segment_entry(owner, i, j, k);
        }else 
            add_segment_entry(owner, i, j, k);
    }
    void layout_table::clean(){
        this->xsegment_count = this->segment_count;
        this->segment_count = 0;
    }
    void layout_table::print(){
        for(int i=0; i < this->segment_count; i++){
            printf("PROFILE: %d; R%d: %d %d\n", this->object->id, this->segment[i].owner, this->segment[i].i, this->segment[i].j);
        }
    }

    void perform_forwarding(p_profile** profiles, size_t count)
    {
        ambient::packets::packet* layout_packet;
        ambient::packets::packet* block_packet;
        for(int k = 0; k < count; k++){
            if(profiles[k]->get_xscope() == NULL || !profiles[k]->get_xscope()->involved()) continue;
            if(ambient::rank() == profiles[k]->get_xmaster()){ // receive and forward new layout
                printf("I'm the old master of the profile (%d)\n", ambient::rank());
                for(int i=0; i < (profiles[k]->get_grid_dim().x*profiles[k]->get_grid_dim().y); i++){
                    layout_packet = ambient::groups::recv<layout_packet_t>(ambient::groups::group_map("ambient"), alloc_t<layout_packet_t>());
                    layout_packet->set(A_DEST_FIELD, profiles[layout_packet->get<int>(A_LAYOUT_P_OP_ID_FIELD)]->
                                                     layout->get_entry(layout_packet->get<int>(A_LAYOUT_P_I_FIELD), 
                                                                       layout_packet->get<int>(A_LAYOUT_P_J_FIELD), 
                                                                       layout_packet->get<int>(A_LAYOUT_P_K_FIELD))->owner); 

                    ambient::groups::send(layout_packet, ambient::groups::group_map("ambient"));
                }
            }
            for(int i=0; i < profiles[k]->layout->xsegment_count; i++){ // receive commands
                layout_packet = ambient::groups::recv<layout_packet_t>(ambient::groups::group_map("ambient"), alloc_t<layout_packet_t>());
                void* header = profiles[layout_packet->get<int>(A_LAYOUT_P_OP_ID_FIELD)]->group(layout_packet->get<int>(A_LAYOUT_P_I_FIELD),
                                                                                                layout_packet->get<int>(A_LAYOUT_P_J_FIELD),
                                                                                                layout_packet->get<int>(A_LAYOUT_P_K_FIELD))->header; 
                block_packet = pack(*profiles[layout_packet->get<int>(A_LAYOUT_P_OP_ID_FIELD)]->packet_type, header, layout_packet->get<int>(A_LAYOUT_P_OWNER_FIELD),
                                                                                                                     "P2P", // communication type
                                                                                                                     k,     // profile id (in terms of profiles array positioning)
                                                                                                                     layout_packet->get<int>(A_LAYOUT_P_I_FIELD),
                                                                                                                     layout_packet->get<int>(A_LAYOUT_P_J_FIELD),
                                                                                                                     layout_packet->get<int>(A_LAYOUT_P_K_FIELD));
                ambient::groups::send(block_packet, ambient::groups::group_map("ambient"));
            }
        }
    }

// exit from this function should be blocking
    void apply_change_set(p_profile** profiles, size_t count)
    {
        ambient::packets::packet* layout_packet;
        ambient::packets::packet* block_packet;
        for(int k = 0; k < count; k++){
            if(ambient::scope.master()){
                for(int i=0; i < (profiles[k]->get_grid_dim().x*profiles[k]->get_grid_dim().y - profiles[k]->layout->segment_count); i++){
                    layout_packet = ambient::groups::recv<layout_packet_t>(ambient::scope.get_group(), alloc_t<layout_packet_t>());
                    profiles[layout_packet->get<int>(A_LAYOUT_P_OP_ID_FIELD)]->layout->update_map_entry(layout_packet->get<int>(A_LAYOUT_P_OWNER_FIELD),
                                                                                                        layout_packet->get<int>(A_LAYOUT_P_I_FIELD)    ,
                                                                                                        layout_packet->get<int>(A_LAYOUT_P_J_FIELD)    ,
                                                                                                        layout_packet->get<int>(A_LAYOUT_P_K_FIELD)    );
                }
            }else{
                for(int i=0; i < profiles[k]->layout->segment_count; i++){
                    layout_packet = pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                          ambient::scope.get_group()->master, 
                                                          "P2P", // communication type
                                                          k,     // profile id (in terms of profiles array positioning
                                                          profiles[k]->layout->segment[i].owner, 
                                                          profiles[k]->layout->segment[i].i, 
                                                          profiles[k]->layout->segment[i].j, 
                                                          profiles[k]->layout->segment[i].k);
                    ambient::groups::send(layout_packet, ambient::scope.get_group());
                }
            }
//            printf("This profile's master is %d (old is %d)\n", profiles[k]->get_master(),profiles[k]->get_xmaster()); 
        } // ok, the layout has been gathered

        for(int k = 0; k < count; k++){
            if(profiles[k]->need_init) profiles[k]->postprocess();
            else{
                if(ambient::scope.master()){ // send the layout to the previous master
                    for(int i=0; i < profiles[k]->get_grid_dim().x*profiles[k]->get_grid_dim().y; i++){
                        layout_packet = pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                              profiles[k]->get_xmaster(), 
                                                              "P2P", // communication type
                                                              k,     // profile id (in terms of profiles array positioning
                                                              profiles[k]->layout->segment[i].owner, 
                                                              profiles[k]->layout->segment[i].i, 
                                                              profiles[k]->layout->segment[i].j,  // not a segment but a map!
                                                              profiles[k]->layout->segment[i].k);
                        ambient::groups::send(layout_packet, ambient::groups::group_map("ambient"));
                    }
                }
                perform_forwarding(profiles, count); // in case if old master/block owners are the part of the selected scope

                for(int i=0; i < profiles[k]->layout->segment_count; i++){
// receive the blocks and add them to profile
                    packet_t* packet_type = profiles[k]->packet_type;
                    block_packet = ambient::groups::recv( *packet_type, ambient::groups::group_map("ambient"), alloc_t(*packet_type) );
                    profiles[block_packet->get<int>(A_BLOCK_P_OP_ID_FIELD)]->group(block_packet->get<int>(A_BLOCK_P_I_FIELD),
                                                                                   block_packet->get<int>(A_BLOCK_P_J_FIELD))->set_memory(block_packet->data);
                }
            }
        }
    }

} }
