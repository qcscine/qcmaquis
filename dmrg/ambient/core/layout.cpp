#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"
#include "ambient/core/auxiliary.h"

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
    layout_table::layout_table(p_profile* profile) 
    : profile(profile), count(0), segment_count(0), request_count(0)
    {
        this->reserved_x = 0;
        this->reserved_y = 0;
        remap();
    }

    void layout_table::remap(){
        int y_size = this->profile->dim.y / (this->profile->get_group_dim().y*this->profile->get_item_dim().y);
        int x_size = this->profile->dim.x / (this->profile->get_group_dim().x*this->profile->get_item_dim().x);
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
    }

    layout_table_entry* layout_table::get_entry(int i, int j, int k){
        if(map[i][j] == NULL) throw race_condition_e(); // to extend for situation when outdated
        return map[i][j];
    }

    layout_table_entry* layout_table::operator()(int i, int j, int k){
        assert(map[i][j] != NULL);
        return map[i][j];
    }

    void layout_table::add_segment_entry(int owner, int i, int j, int k){
        if(segment_count == segment.size()) segment.resize(segment_count+1);
        segment[segment_count++] = layout_table_entry(owner, i, j, k);
    }

    void layout_table::add_request_entry(int i, int j, int k){
        if(request_count == requests.size()) requests.resize(request_count+1);
        requests[request_count++] = layout_table_entry(-1, i, j, k);
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
    void layout_table::record(int i, int j, int k) {
        for(int s=0; s < this->segment_count; s++)
            if(this->segment[s].i == i &&
               this->segment[s].j == j &&
               this->segment[s].k == k) return; // avoiding redunant information // that is - hangs in mpi

        this->profile->group(i,j,k)->owner = ambient::rank();
        add_segment_entry(ambient::rank(), i, j, k);
    }
    void layout_table::request(int i, int j, int k){
        if(this->profile->state == COMPOSING) 
            return record(i,j,k);
        for(int s=0; s < this->request_count; s++)
            if(this->requests[s].i == i &&
               this->requests[s].j == j &&
               this->requests[s].k == k) return; // avoiding redunant requests
        add_request_entry(i, j, k);
    }
    void layout_table::clean(){
        this->segment_count = 0;
        this->request_count = 0;
    }
    void layout_table::print(){
        for(int i=0; i < this->segment_count; i++){
            printf("PROFILE: %d; R%d: %d %d\n", this->profile->id, this->segment[i].owner, this->segment[i].i, this->segment[i].j);
        }
    }

    ambient::packets::packet* package(p_profile* profile, const char* state, int i, int j, int k, int dest)
    {
        if(*state == 'P'){
             if(profile->associated_proxy == NULL) throw race_condition_e();
             profile = profile->associated_proxy; // GLOBAL REDUCTION HANDLING
        }
        void* header = profile->group(i, j, k)->header; 
        if(header == NULL) throw race_condition_e(); // to extend for situation when outdated
        return pack(*profile->packet_type, header, dest, "P2P", *profile->group_id, profile->id, state, i, j, k, NULL);
    }

    void forward_block(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->profile;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'R') return; // REQUEST TRANSFER TO THE NEW OWNER ACTION
        if(!profile->xinvolved()) return;
        try{
            in_q.manager->emit(package(profile, (const char*)pack->get(A_LAYOUT_P_STATE_FIELD), pack->get<int>(A_LAYOUT_P_I_FIELD), pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                       pack->get<int>(A_LAYOUT_P_K_FIELD), pack->get<int>(A_LAYOUT_P_OWNER_FIELD)));
        }catch(race_condition_e){
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    void forward_layout(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->profile;
        if(!profile->xinvolved()) return; // can be omitted I guess
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM X OWNER ACTION
        try{
            layout_table_entry* entry = profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD), 
                                                                   pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                                                   pack->get<int>(A_LAYOUT_P_K_FIELD));
            pack->set(A_DEST_FIELD, entry->get_owner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
            in_q.manager->emit(pack);
        }catch(race_condition_e){
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    void update_layout(packet_manager::typed_q& in_q)
    {
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->profile;
        if(!profile->get_scope()->is_master()) return;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'U' &&       // UPDATE ACTION
           pack->get<char>(A_LAYOUT_P_ACTION) != 'C' ) return;
        profile->layout->update_map_entry(pack->get<int>(A_LAYOUT_P_OWNER_FIELD),
                                          pack->get<int>(A_LAYOUT_P_I_FIELD)    ,
                                          pack->get<int>(A_LAYOUT_P_J_FIELD)    ,
                                          pack->get<int>(A_LAYOUT_P_K_FIELD)    );
        if(pack->get<char>(A_LAYOUT_P_ACTION) == 'C') return; // COMPOSE ACTION

        if(profile->get_xmaster() == profile->get_master()){
            pack->set(A_DEST_FIELD, profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                               pack->get<int>(A_LAYOUT_P_J_FIELD))->get_xowner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
        }else{ 
            pack->set(A_DEST_FIELD, profile->get_xmaster());
            pack->set(A_LAYOUT_P_ACTION, "INFORM X OWNER");
        }
        in_q.manager->emit(pack);
    }

    void apply_changes(p_profile** profiles, size_t count)
    {
        for(int k = 0; k < count; k++){
            const char* action = "UPDATE";
            if(profiles[k]->state == COMPOSING){
                profiles[k]->postprocess(); 
                action = "COMPOSE";
            }
            for(int i=0; i < profiles[k]->layout->segment_count; i++){
                world()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                                   scope.get_master_g(), "P2P", action,
                                                                  *profiles[k]->group_id, profiles[k]->id, "GENERIC",
                                                                   profiles[k]->layout->segment[i].owner, 
                                                                   profiles[k]->layout->segment[i].i, 
                                                                   profiles[k]->layout->segment[i].j, 
                                                                   profiles[k]->layout->segment[i].k));
            }
            for(int i=0; i < profiles[k]->layout->request_count; i++){
                world()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                                   profiles[k]->get_master(), "P2P", 
                                                                  "INFORM OWNER ABOUT REQUEST",
                                                                  *profiles[k]->group_id, profiles[k]->id, "GENERIC",
                                                                   ambient::rank(), // forward target
                                                                   profiles[k]->layout->requests[i].i, 
                                                                   profiles[k]->layout->requests[i].j, 
                                                                   profiles[k]->layout->requests[i].k));
            }
        }
    }

} }
