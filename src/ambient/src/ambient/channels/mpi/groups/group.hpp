#include "ambient/utils/ceil.h"
#include "ambient/utils/hashmap.hpp"

// goals: 
// 1. minimize the number of transfers
// 2. unleash most of the parallelism *
// 3. balance memory consumption
// 4. balance comp. resources usage *

// methods to reduce size of id:
// 1. discrete regions of procs
// 2. minus mask
// 3. ranges of procs (supplementary array?)

// bet factors:
// 1. dependency number
// 2. size of the matrix to transfer
// 3. cost of the operation

// GOAL:
// MIN(TRANSFER TIME) + MIN(COMPUTATIONAL TIME)

// MIN(TRANSFER TIME):
// OBJECTS ARE STICKY (MULTIPLE READ COPIES?)
// CAREFULL CHOOSING OF OPEATIONS LOCATIONS

/*

GRAPH PROBLEM:

1 2 3 4 ... N

{ 1 .. n11 }
{ 1 .. n12 }
{ 1 .. n13 }
{ 1 .. n1m }

{ 2 .. n21 }
{ 2 .. n22 }
{ 2 .. n2k }

...


.  .  .  .  .  .  .
.  .  .  .  .  .  .
.  .  .  .  .  .  .
.  .  .  .  .  .  .
.  .  .  .  .  .  .
.  .  .  .  .  .  .
.  .  .  .  .  .  .


I. resource graph consumption:
for each round:

1. estimate total compute time / number of procs // move ctxt_select to templates or push ?
2. sort operations according to the number of rounds (dependencies) or other metrics
3. find nearest available resource entry and eat it


II. operation graph consumption:
for each round:

1. estimate nearest entry point
2. estimate nearest next point
3. eat up the available amount of time
*/

namespace ambient { namespace channels { namespace mpi {

    inline group::~group(){
        free(this->vacations);
        free(this->id.first);
        free(this->members);
        free(this->members_g);
    }

    inline group::group(const char* name, int master, MPI_Comm parent)
    : members(NULL), vacant_level(0), occupancy(false)
    {
        this->parent = NULL;
        this->mpi_comm = parent;
        MPI_Comm_group(this->mpi_comm, &this->mpi_group);
        MPI_Group_size(this->mpi_group, &this->count);
        MPI_Group_rank(this->mpi_group, &this->rank);
        this->vacations = (int*)malloc(sizeof(int)*this->count);
        memset(this->vacations, 0, sizeof(int)*this->count);
        this->members_g = (int*)malloc(sizeof(int)*this->count);
        for(int i=0; i < this->count; i++) this->members_g[i] = i;
        if(name == NULL) this->name = (char*)"tmp";
        else this->name = name;
        this->master = master;
        this->id = hash_gid();
        group_id_map(this->id, this);
    }

    inline std::pair<size_t*,size_t> group::hash_gid(){ // now we don't really need this
        int index;
        int old_id_len;
        int id_len = 0;
        int rank_g;
        int volume;
        size_t* hash_id = NULL;

        MPI_Comm_size(MPI_COMM_WORLD, &volume);
        if(this->count == volume){
            hash_id = (size_t*)malloc(sizeof(size_t));
            id_len = *hash_id = 1; // 1 as a first bit stands for mirroring of all values
        }else{
            for(int i = 0; i < this->count; i++){
                rank_g = this->translate_up_rank(i)+1;
                if(rank_g >= id_len*32){
                    old_id_len = id_len;
                    id_len = __a_ceil(rank_g/32);
                    hash_id = (size_t*)realloc(hash_id, sizeof(size_t)*id_len);
                    while(old_id_len < id_len) hash_id[old_id_len++] = 0; // memset with 0
                }
                index = rank_g/32;
                hash_id[index] |= 1 << rank_g % 32;
            }
        }
        return std::pair<size_t*,size_t>(hash_id,1);
    }

    inline group::group(const char* name, int master, group* parent)
    : count(0), members(NULL), members_g(NULL), vacations(NULL), vacant_level(0), occupancy(false)
    {
        this->parent = parent;
        if(name == NULL) this->name = (char*)"tmp";
        else this->name = name;
        this->master = master;
        this->parent->children.insert(this);
    }

    inline void group::add(const int* procs, int count, bool excl){
        if(count <= 0) return;
        this->members   = (int*)realloc(this->members,   (this->count+count)*sizeof(int));
        this->members_g = (int*)realloc(this->members_g, (this->count+count)*sizeof(int));
        memcpy(&(this->members[this->count]), procs, count*sizeof(int));
        for(int i=0; i < count; i++)
            this->members_g[this->count+i] = this->parent->translate_up_rank(this->members[this->count+i]);
        if(excl == true){
            // TODO: discard existing groups
        }
        this->count += count;
    }

    inline void group::add(int count, bool excl){
        if(count <= 0) return;
        count = std::min(count, parent->count); // avoiding over-allocation problems
        this->members   = (int*)realloc(this->members,   (this->count+count)*sizeof(int));
        this->members_g = (int*)realloc(this->members_g, (this->count+count)*sizeof(int));
        for(int i=0; i < count; i++){ 
            this->members  [this->count+i] = this->parent->get_vacant();
            this->members_g[this->count+i] = this->parent->translate_up_rank(this->members[this->count+i]);
        }
        if(excl == true){
            // TODO: discard existing groups
        }
        this->count += count;
    }

    inline int group::get_vacant(){
        for(int i=0; i < this->count; i++){
            if(this->vacations[i] == this->vacant_level){
                this->vacations[i]++;
                return i;
            }
        }
        this->vacant_level++;
        return this->get_vacant(); 
// note:recursion will be endless in case of 0-group
    }

    inline void group::add_range(int first, int last, bool excl){
        int count = last - first + 1;
        int procs[count];
        for(int i=0; i<count; i++) procs[i] = first + i;
        this->add(procs, count, excl);
    }

    inline void group::add_every(int nth, bool excl){
        int count = (int)(this->parent->count / nth);
        int procs[count];
        if(this->parent == NULL) printf("Warning: attempting to modify root group.\n");
        for(int i=0; i < count; i++) procs[i] = (i+1)*nth-1;
        this->add(procs, count, excl);
    }

    inline void group::add_every(bool(*include)(int k), bool excl){
        int procs[this->parent->count];
        int count = 0;
        for(int i=0; i<this->parent->count; i++){
            if(include(i) == true) procs[count++] = i;
        }
        this->add(procs, count, excl);
    }

    inline void group::add_every(bool(*include)(int k), bool(*excl)(int k)){
        int procs[this->parent->count];
        int count = 0;
        bool exclusive = false;

        for(int i=0; i<this->parent->count; i++){
            if(include(i)){
                if(exclusive != excl(i)){
                    if(count) this->add(procs, count, exclusive);
                    exclusive = !exclusive;
                    count = 0;
                }
                procs[count++] = i;
            }
        }
    }

    inline void group::add_intersection(const group* b, int* count){
        for(int i=0; i < this->parent->count; i++){
            for(int j=0; j < b->count; j++){
                if(this->parent->members_g[i] == b->members_g[j]){
                    if(count != NULL && (*count)-- <= 0) return;
                    this->add(&i, 1);
                    break;
                }
            }
        }
    }

    inline void group::add_every_intersection(const group* b, int nth, int* count){
        assert(false);
        int countdown = nth;
        for(int i=0; i < this->parent->count; i++){
            for(int j=0; j < b->count; j++){
                if(this->parent->members_g[i] == b->members_g[j]){
                    if(count != NULL && (*count)-- <= 0) return;
                    if(--countdown == 0){ this->add(&i, 1); countdown = nth; }
                    break;
                }
            }
        }
    }

    inline void group::add_substraction(const group* b, int* count){
        for(int i=0; i < this->parent->count; i++){
            for(int j=0; j < b->count; j++){
                if(this->parent->members_g[i] == b->members_g[j]) break;
                if(j == b->count-1){
                    if(count != NULL && (*count)-- <= 0) return;
                    this->add(&i, 1);
                }
            }
        }
    }

    inline void group::reorder(int* new_ranks){ //  rank i becomes new_ranks[i];
        int* members_new = (int*)malloc(this->count*sizeof(int));
        for(int i=0; i < this->count; i++){
            members_new[new_ranks[i]] = this->members[i];
        }
        free(this->members);
        this->members = members_new;

        for(std::set<group*>::iterator it = this->children.begin(); it != this->children.end(); it++ )
            for(int i=0; i < (*it)->count; i++)
                (*it)->members[i] = new_ranks[(*it)->members[i]];
    }

    inline void group::reorder(int(*order)(int r)){
        int* members_new = (int*)malloc(this->count*sizeof(int));
        for(int i=0; i < this->count; i++){
            members_new[order(i)] = this->members[i];
        }
        free(this->members);
        this->members = members_new;

        for(std::set<group*>::iterator it = this->children.begin(); it != this->children.end(); it++ )
            for(int i=0; i < (*it)->count; i++)
                (*it)->members[i] = order((*it)->members[i]);
    }

// note: be sure to feed this with the rank inside group (i) not members[i]
// as members[i] contains translation to the parent group
    inline int group::translate_up_rank(int rank, group* parent) const {
        int rank_n = rank;
        const group* iterator = this;
        if(!parent) parent = group_map("ambient");
        assert(rank_n >= 0);
        if(parent == this) return rank_n;
        do{
            if(rank_n >= iterator->count || (iterator->parent == NULL && parent != NULL)){
                printf("Warning: the rank doesn't belongs to this group\n"); 
                break; 
            }
            rank_n = iterator->members[rank_n];
            iterator = iterator->parent;
        }while(iterator != parent);
        assert(rank_n >= 0);
        return rank_n;
    }

// note: be sure to feed this with the rank inside group (i) not members[i]
// as members[i] contains translation to the parent group
// note: invalid reads unsafe...
    inline int group::translate_down_rank(int rank, group* child) const {
        int rank_n = rank;
        const group* iterator = child;
        if(child == this) return rank;
        else{
            rank_n = group::translate_down_rank(rank, child->parent);
            for(int i=0; i < child->count; i++)
            if(child->members[i] == rank_n) return i;
        }
        assert(false);
        return -1;
    }

    inline void group::commit(){
        this->id = hash_gid();
        group* grp = group_id_map(this->id);
        if(grp != NULL) throw grp;
        else group_id_map(this->id, this);

        MPI_Group_incl(this->parent->mpi_group, this->count, this->members, &this->mpi_group);
        MPI_Group_rank(this->mpi_group, &this->rank);
        this->vacations = (int*)malloc(sizeof(int)*this->count);
        memset(this->vacations, 0, sizeof(int)*this->count);
    }

    inline int group::get_master(){
        return this->master;
    }

    inline int group::get_master_g(){
        return translate_up_rank(this->master);
    }

    inline int group::get_rank(){
        return this->rank;
    }

    inline bool group::involved(){
        return (this->rank != MPI_UNDEFINED);
    }

    inline bool group::is_master(){
        return (this->rank == this->master);
    }

    inline bool group::occupied(){
        return this->occupancy;
    }

    inline void group::occupy(){
        this->occupancy = true;
    }

    inline void group::idle(){
        this->occupancy = false; 
    }

    inline size_t group::get_size(){
        return this->count;
    }

    inline const char* group::get_name(){
        return this->name;
    }

    inline size_t group::get_member(size_t i){
        assert(i < this->get_size());
        return this->members_g[i];
    }

    inline group* group_map(const char* name, group* instance){ // timings: 0.36 seconds
        static std::map<std::string,group*> map;
        if(instance != NULL){
            if(map.find(name) != map.end()) map.find(name)->second = instance; // can delete the group here
            else map.insert(std::pair<std::string,group*>(name,instance));
            return instance;
        }
        if(map.find(name) == map.end()) instance = NULL; // wasn't able to find requested group
        else instance = map.find(name)->second;
        return instance;
    }

    inline group* group_id_map(std::pair<size_t*,size_t> id, group* instance){
        static fhashmap map;
        group** original = (group**)map.get(id.first, id.second);
        if(instance != NULL){ 
            *original = instance; 
            group_map(instance->name, instance); 
        }
        return *original;
    }

} } }
