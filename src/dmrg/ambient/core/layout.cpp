#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"

#include "ambient/core/operation/operation.h"

namespace ambient{
#include "ambient/core/operation/operation.pp.sa.hpp"
}

namespace ambient{ namespace core{

    using namespace ambient::groups;

    layout_table_entry::layout_table_entry()
    { }

    layout_table_entry::layout_table_entry(int owner, int i, int j, int k)
    : xowner(-1), owner(owner), i(i), j(j), k(k) { }

    int layout_table_entry::get_xowner(){
        return (this->xowner == -1 ? this->owner : this->xowner);
    }
    int layout_table_entry::get_owner(){
        return this->owner;
    }

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
        if(map[i][j] == NULL){ 
            this->count++;
            map[i][j] = new layout_table_entry(owner, i, j, k);
        }else{
            map[i][j]->xowner = map[i][j]->owner;
            map[i][j]->owner = owner;
        }
    }
    void layout_table::record(int owner, int i, int j, int k) {
        for(int s=0; s < this->segment_count; s++)
            if(this->segment[s].i == i &&
               this->segment[s].j == j &&
               this->segment[s].k == k) return; // avoiding redunant information // that is - hangs in mpi

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


// REDISTRIBUTION HEAVY STUFF //

    ambient::packets::packet* package(p_profile* profile, int id, int i, int j, int k, int dest){
        void* header = profile->group(i, j, k)->header; 
        return pack(*profile->packet_type, header, dest, "P2P", id, i, j, k, NULL);
    }

    void forward_layout(p_profile**& profiles, packet_manager::typed_q& in_q, packet_manager::typed_q& out_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        layout_table_entry* entry = profiles[pack->get<int>(A_LAYOUT_P_OP_ID_FIELD)]->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD), 
                                                                                                        pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                                                                                        pack->get<int>(A_LAYOUT_P_K_FIELD));
        pack->set(A_DEST_FIELD, entry->get_owner());
        out_q.push(pack);
    }

    void forward_block(p_profile**& profiles, packet_manager::typed_q& in_q, packet_manager::typed_q& out_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        out_q.push(package(profiles[pack->get<int>(A_LAYOUT_P_OP_ID_FIELD)], pack->get<int>(A_LAYOUT_P_OP_ID_FIELD), 
                           pack->get<int>(A_LAYOUT_P_I_FIELD), pack->get<int>(A_LAYOUT_P_J_FIELD), pack->get<int>(A_LAYOUT_P_K_FIELD), 
                           pack->get<int>(A_LAYOUT_P_OWNER_FIELD)));
    }

    void update_layout(p_profile**& profiles, packet_manager::typed_q& in_q, packet_manager::typed_q& out_q)
    {
        ambient::packets::packet* pack = in_q.get_target_packet();
        int k = pack->get<int>(A_LAYOUT_P_OP_ID_FIELD);
        profiles[k]->layout->update_map_entry(pack->get<int>(A_LAYOUT_P_OWNER_FIELD),
                                             pack->get<int>(A_LAYOUT_P_I_FIELD)    ,
                                             pack->get<int>(A_LAYOUT_P_J_FIELD)    ,
                                             pack->get<int>(A_LAYOUT_P_K_FIELD)    );
        if(profiles[k]->need_init) return;

        if(profiles[k]->get_xmaster() == profiles[k]->get_master())
            pack->set(A_DEST_FIELD, profiles[k]->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                                   pack->get<int>(A_LAYOUT_P_J_FIELD))->get_xowner());
        else 
            pack->set(A_DEST_FIELD, profiles[k]->get_xmaster());
        out_q.push(pack);
    }

    void integrate_block(p_profile**& profiles, packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        profiles[pack->get<int>(A_BLOCK_P_OP_ID_FIELD)]->group(pack->get<int>(A_BLOCK_P_I_FIELD), 
                                                               pack->get<int>(A_BLOCK_P_J_FIELD))->set_memory(pack->data);
    }

    void perform_forwarding(p_profile** profiles, size_t count)
    {
        packet_manager::typed_q* layout_in_q  = world()->get_manager()->add_typed_q(get_t<layout_packet_t>(), packet_manager::IN, 30);
        packet_manager::typed_q* layout_out_q = world()->get_manager()->add_typed_q(get_t<layout_packet_t>(), packet_manager::OUT, 30);
        packet_manager::typed_q* block_out_q[count];

        for(int k = 0; k < count; k++){
            if(!profiles[k]->xinvolved()) continue;
            if(ambient::rank() == profiles[k]->get_xmaster()){
                layout_in_q->packet_delivered += new core::operation(forward_layout, &profiles, layout_in_q, layout_out_q);
                break;
            }
        }
        for(int k = 0; k < count; k++){
            if(!profiles[k]->xinvolved()) continue;
            block_out_q[k] = world()->get_manager()->add_typed_q(*profiles[k]->packet_type, packet_manager::OUT, 30);
            layout_in_q->packet_delivered += new core::operation(forward_block, &profiles, layout_in_q, block_out_q[k]);
            break;
        }
    }

    void apply_change_set(p_profile** profiles, size_t count)
    {
        packet_manager::typed_q* layout_in_q;
        packet_manager::typed_q* layout_out_q;
        packet_manager::typed_q* block_in_q[count]; 

        layout_out_q = world()->get_manager()->add_typed_q(get_t<layout_packet_t>(), packet_manager::OUT, 30);
        if(scope.master()){
            layout_in_q  = world()->get_manager()->add_typed_q(get_t<layout_packet_t>(), packet_manager::IN, 30);
            layout_in_q->packet_delivered += new core::operation(update_layout, &profiles, layout_in_q, layout_out_q);
        }

        for(int k = 0; k < count; k++)
            for(int i=0; i < profiles[k]->layout->segment_count; i++)
                layout_out_q->push(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                         scope.get_group()->get_master_g(), 
                                                         "P2P", // communication type: peer to peer
                                                         k,     // profile id in terms of profiles-array pos
                                                         profiles[k]->layout->segment[i].owner, 
                                                         profiles[k]->layout->segment[i].i, 
                                                         profiles[k]->layout->segment[i].j, 
                                                         profiles[k]->layout->segment[i].k));
        for(int k = 0; k < count; k++)
            if(profiles[k]->need_init) 
                profiles[k]->postprocess();
            else{
                block_in_q[k] = world()->get_manager()->add_typed_q(*profiles[k]->packet_type, packet_manager::IN, 30);
                block_in_q[k]->packet_delivered += new core::operation(integrate_block, &profiles, block_in_q[k]);
            }
        perform_forwarding(profiles, count);
        world()->get_manager()->process();
    }

} }
