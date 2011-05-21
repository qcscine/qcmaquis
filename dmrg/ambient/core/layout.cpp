#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"
#include "ambient/core/auxiliary.h"

namespace ambient{ namespace core{

    using namespace ambient::groups;

    layout_table::composite_marker::composite_marker(){
        this->active = true;
        this->imarker = 0;
        this->jmarker = 0;
    }
    void layout_table::composite_marker::clear(){
        this->active = false;
    }
    void layout_table::composite_marker::mark(int i, int j){
        this->active = true;
        this->imarker = i;
        this->jmarker = j;
    }
    bool layout_table::composite_marker::has_marked(int i, int j){
        if(!this->active) return false;
        if(i >= this->imarker || j >= this->jmarker) return true;
        return false;
    }

    layout_table::entry::entry()
    { }

    layout_table::entry::entry(int owner, int i, int j)
    : xowner(-1), owner(owner), i(i), j(j) { }

    int layout_table::entry::get_xowner(){
        return (this->xowner == -1 ? this->owner : this->xowner);
    }
    int layout_table::entry::get_owner(){
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
        int y_size = __a_ceil(this->profile->get_dim().y / this->profile->get_mem_t_dim().y);
        int x_size = __a_ceil(this->profile->get_dim().x / this->profile->get_mem_t_dim().x);
        if(this->reserved_x >= x_size && this->reserved_y >= y_size) return;
        for(int i = 0; i < y_size; i++){
            if(i >= this->reserved_y) map.push_back(std::vector<entry*>());
            for(int j = 0; j < x_size; j++){
                if(j >= this->reserved_x || i >= this->reserved_y) 
                    map[i].push_back(NULL); //new entry(-1, i, j);
            }
        }
        if(x_size > this->reserved_x) this->reserved_x = x_size;
        if(y_size > this->reserved_y) this->reserved_y = y_size;
    }

    std::vector<core::layout_table::entry>& layout_table::get_list(){
        return (this->segment_count != 0 ? this->segment : this->requests);
    }

    layout_table::entry* layout_table::get_entry(int i, int j){
        if(map[i][j] == NULL) throw race_condition_e(); // to extend for situation when outdated
        return map[i][j];
    }

    layout_table::entry* layout_table::operator()(int i, int j){
        assert(map[i][j] != NULL);
        return map[i][j];
    }

    void layout_table::add_segment_entry(int owner, int i, int j){
        if(segment_count == segment.size()) segment.resize(segment_count+1);
        segment[segment_count++] = entry(owner, i, j);
    }

    void layout_table::add_request_entry(int i, int j){
        if(request_count == requests.size()) requests.resize(request_count+1);
        requests[request_count++] = entry(-1, i, j);
    }

    void layout_table::update_map_entry(int owner, int i, int j){
        if(map[i][j] == NULL){ 
            this->count++;
            map[i][j] = new entry(owner, i, j);
        }else{
            map[i][j]->xowner = map[i][j]->owner;
            map[i][j]->owner = owner;
        }
    }
    void layout_table::record(int i, int j) {
        for(int s=0; s < this->segment_count; s++)
            if(this->segment[s].i == i &&
               this->segment[s].j == j) return; // avoiding redunant information // that is - hangs in mpi

        this->profile->block(i,j)->owner = ambient::rank();
        add_segment_entry(ambient::rank(), i,j);
    }
    void layout_table::request(int i, int j){
        if(this->init_marker.has_marked(i,j)){
            return record(i,j);
        }
        for(int s=0; s < this->request_count; s++)
            if(this->requests[s].i == i &&
               this->requests[s].j == j) return; // avoiding redunant requests
        add_request_entry(i, j);
    }
    void layout_table::clean(){
        this->segment_count = 0;
        this->request_count = 0;
        this->segment.resize(0);
        this->requests.resize(0);
    }
    void layout_table::print(){
        for(int i=0; i < this->segment_count; i++){
            printf("PROFILE: %d; R%d: %d %d\n", this->profile->id, this->segment[i].owner, this->segment[i].i, this->segment[i].j);
        }
    }

    ambient::packets::packet* package(p_profile* profile, const char* state, int i, int j, int dest)
    {
        if(*state == 'P'){
             if(profile->associated_proxy == NULL) throw race_condition_e();
             profile = profile->associated_proxy; // GLOBAL REDUCTION HANDLING
        }
        void* header = profile->block(i,j)->header;
        if(header == NULL){ printf("Failing on id %d, getting %d %d (dim: %d x %d)\n", (int)profile->id, i, j, profile->dim.y, profile->dim.x); assert(false); printf("Warning: no header available\n"); throw race_condition_e(); }// to extend for situation when outdated
        return pack(*profile->packet_type, header, dest, "P2P", *profile->group_id, profile->id, state, i, j, NULL);
    }

    void forward_block(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->profile;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'R') return; // REQUEST TRANSFER TO THE NEW OWNER ACTION
        if(!profile->xinvolved()) return;
        try{
            in_q.manager->emit(package(profile, (const char*)pack->get(A_LAYOUT_P_STATE_FIELD), 
                                       pack->get<int>(A_LAYOUT_P_I_FIELD), pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                       pack->get<int>(A_LAYOUT_P_OWNER_FIELD)));
        }catch(race_condition_e){
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    void forward_layout(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->profile;
        if(!profile->xinvolved()) return; // can be omitted I guess
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM X OWNER ACTION
        try{
            layout_table::entry* entry = profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD), 
                                                                    pack->get<int>(A_LAYOUT_P_J_FIELD)); 
            pack->set(A_DEST_FIELD, entry->get_owner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
            in_q.manager->emit(pack);
        }catch(race_condition_e){
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
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
                                          pack->get<int>(A_LAYOUT_P_J_FIELD)    );
        if(pack->get<char>(A_LAYOUT_P_ACTION) == 'C') return; // COMPOSE ACTION

        if(profile->get_xmaster() == profile->get_master()){
            pack->set(A_DEST_FIELD, profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                               pack->get<int>(A_LAYOUT_P_J_FIELD))->get_xowner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
        }else{
            pack->set(A_DEST_FIELD, profile->get_xmaster());
            pack->set(A_LAYOUT_P_ACTION, "INFORM X OWNER");
            assert(pack->get<int>(A_DEST_FIELD) >= 0);
        }
        in_q.manager->emit(pack);
    }

    void apply_changes(p_profile** profiles, size_t count)
    {
        for(int k = 0; k < count; k++){
            for(int i=0; i < profiles[k]->layout->segment_count; i++){
                const char* action = "UPDATE";
                if(profiles[k]->layout->init_marker.active){
                    if(profiles[k]->layout->init_marker.
                                            has_marked(profiles[k]->layout->segment[i].i, 
                                                       profiles[k]->layout->segment[i].j))
                    {
                        action = "COMPOSE";
                        int ii = profiles[k]->layout->segment[i].i;
                        int jj = profiles[k]->layout->segment[i].j;
                        profiles[k]->postprocess(ii,jj);
                    }
                }
                world()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                                   scope.get_master_g(), "P2P", action,
                                                                  *profiles[k]->group_id, profiles[k]->id, "GENERIC",
                                                                   profiles[k]->layout->segment[i].owner, 
                                                                   profiles[k]->layout->segment[i].i, 
                                                                   profiles[k]->layout->segment[i].j));
            }
            for(int i=0; i < profiles[k]->layout->request_count; i++){
                if(profiles[k]->block(profiles[k]->layout->requests[i].i, 
                                      profiles[k]->layout->requests[i].j)->available()) continue; // avoiding redunant requests
                world()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                                   profiles[k]->get_master(), "P2P", 
                                                                  "INFORM OWNER ABOUT REQUEST",
                                                                  *profiles[k]->group_id, profiles[k]->id, "GENERIC",
                                                                   ambient::rank(), // forward target
                                                                   profiles[k]->layout->requests[i].i, 
                                                                   profiles[k]->layout->requests[i].j));
            }
        }
    }

} }
