#include "ambient/core/layout.h"
#include "ambient/ambient.h"
#include "ambient/core/p_profile.h"
#include "ambient/packets/auxiliary.hpp"
#include "ambient/groups/auxiliary.hpp"
#include "ambient/groups/packet_manager.h"
#include "ambient/core/operation/operation.h"
#include "ambient/core/operation/operation.pp.sa.hpp"

namespace ambient{ namespace core{

    using namespace ambient::groups;

    class race_condition_e{ };

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

        if(scope.is_master()){ 
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

    ambient::packets::packet* package(p_profile* profile, int i, int j, int k, int dest){
        void* header = profile->group(i, j, k)->header; 
        assert(header != NULL);
        return pack(*profile->packet_type, header, dest, "P2P", *profile->group_id, profile->id, i, j, k, NULL);
    }

    void forward_block(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->object;
        if(!profile->xinvolved()) return;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'R') return; // REQUEST TRANSFER TO THE NEW OWNER ACTION
        in_q.manager->emit(package(profile, pack->get<int>(A_LAYOUT_P_I_FIELD), pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                   pack->get<int>(A_LAYOUT_P_K_FIELD), pack->get<int>(A_LAYOUT_P_OWNER_FIELD)));
    }

    void forward_layout(packet_manager::typed_q& in_q){
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->object;
        if(!profile->xinvolved()) return; // can be omitted I guess
        if(ambient::rank() != profile->get_xmaster()) return;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'I') return; // INFORM X OWNER ACTION

        try{
            layout_table_entry* entry = profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD), 
                                                                   pack->get<int>(A_LAYOUT_P_J_FIELD), 
                                                                   pack->get<int>(A_LAYOUT_P_K_FIELD));
            pack->set(A_DEST_FIELD, entry->get_owner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
            in_q.manager->emit(pack);
            printf("R%d requeting piece of layout %d %d %d of %u:%d from %d\n", ambient::rank(),pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                             pack->get<int>(A_LAYOUT_P_J_FIELD),
                                                             pack->get<int>(A_LAYOUT_P_K_FIELD),
	    						     *(unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD),
							     pack->get<int>(A_LAYOUT_P_ID_FIELD), pack->get<int>(A_DEST_FIELD));

        }catch(race_condition_e){
            in_q.manager->emit(pack); // re-throwing the packet for future handling
        }
    }

    void update_layout(packet_manager::typed_q& in_q)
    {
        ambient::packets::packet* pack = in_q.get_target_packet();
        p_profile* profile = p_profile_map.find((unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD), 1, pack->get<int>(A_LAYOUT_P_ID_FIELD))->object;
        if(!profile->get_scope()->is_master()) return;
        if(pack->get<char>(A_LAYOUT_P_ACTION) != 'U' &&       // UPDATE ACTION
           pack->get<char>(A_LAYOUT_P_ACTION) != 'C' ) return;
        profile->layout->update_map_entry(pack->get<int>(A_LAYOUT_P_OWNER_FIELD),
                                          pack->get<int>(A_LAYOUT_P_I_FIELD)    ,
                                          pack->get<int>(A_LAYOUT_P_J_FIELD)    ,
                                          pack->get<int>(A_LAYOUT_P_K_FIELD)    );
        if(pack->get<char>(A_LAYOUT_P_ACTION) == 'C'){ printf("Actually I did caught that for profile %d\n", profile->id); return; } // COMPOSE ACTION

        if(profile->get_xmaster() == profile->get_master()){
            pack->set(A_DEST_FIELD, profile->layout->get_entry(pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                               pack->get<int>(A_LAYOUT_P_J_FIELD))->get_xowner());
            pack->set(A_LAYOUT_P_ACTION, "REQUEST TRANSFER TO THE NEW OWNER");
        }else{ 
            pack->set(A_DEST_FIELD, profile->get_xmaster());
            pack->set(A_LAYOUT_P_ACTION, "INFORM X OWNER");
            printf("R%d informing old owner of layout %d %d %d of %u:%d - now %d\n", ambient::rank(),pack->get<int>(A_LAYOUT_P_I_FIELD),
                                                             pack->get<int>(A_LAYOUT_P_J_FIELD),
                                                             pack->get<int>(A_LAYOUT_P_K_FIELD),
	    						     *(unsigned int*)pack->get(A_LAYOUT_P_GID_FIELD),
							     pack->get<int>(A_LAYOUT_P_ID_FIELD), pack->get<int>(A_LAYOUT_P_OWNER_FIELD));
        }
        in_q.manager->emit(pack);
    }

    void apply_changes(p_profile** profiles, size_t count)
    {
        for(int k = 0; k < count; k++){
            const char* action = "UPDATE";
            if(profiles[k]->need_init){
                profiles[k]->postprocess(); 
                action = "COMPOSE";
            }
            for(int i=0; i < profiles[k]->layout->segment_count; i++){
                world()->get_manager()->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), 
                                                                   scope.get_master_g(), "P2P", action,
                                                                  *profiles[k]->group_id, profiles[k]->id,
                                                                   profiles[k]->layout->segment[i].owner, 
                                                                   profiles[k]->layout->segment[i].i, 
                                                                   profiles[k]->layout->segment[i].j, 
                                                                   profiles[k]->layout->segment[i].k));
            }
        }
    }

} }
