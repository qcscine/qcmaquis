#include "ambient/ambient.h"
#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"

namespace ambient{ namespace groups {

    comm_map::comm_map& comm_map::instance()
    {
        static comm_map* singleton = NULL;
        if(!singleton) singleton = new comm_map();
        return *singleton;
    }
    comm_map::comm_map():content(HASH_MAP_PARTITION_SIZE){ }
    
    group** comm_map::get(unsigned int* hash, unsigned int hash_len, int shift) const
    {
        unsigned int hash_w = hash[0] >> shift;
        if(hash_w >= HASH_MAP_PARTITION_SIZE){
            unsigned int hash_cut = (unsigned int)(unsigned char)hash_w; // first log2 of HASH_MAP_PARTITION_SIZE bits
            if(this->content[hash_cut].first == NULL){
                this->content[hash_cut].first = new comm_map();
            }
            return this->content[hash_cut].first->get(hash, hash_len, shift+HASH_MAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            if(this->content[hash_w].first == NULL){
                this->content[hash_w].first = new comm_map();
            }
            return this->content[hash_w].first->get(&hash[1], hash_len-1); 
        }
        return &this->content[hash_w].second;
    }

    group::group(const char* name, int master, MPI_Comm parent): members(NULL), object_count(0), vacant_level(0)
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
        this->name = name;
        this->master = master;
        this->manager = new packet_manager(this);
        this->id = hash_group_id();
        group** original = mpi_comm_map.get(this->id.first, this->id.second);
       *original = this;
        group_map(this->name, this);
    }

// methods to reduce size of id:
// 1. discrete regions of procs
// 2. minus mask
// 3. ranges of procs (supplementary array?)
    std::pair<unsigned int*,size_t> group::hash_group_id()
    {
        int index;
        int old_id_len;
        int id_len = 0;
        int rank_g;
        unsigned int* hash_id = NULL;

        if(this->count == ambient::size()){
            hash_id = (unsigned int*)malloc(sizeof(unsigned int));
            id_len = *hash_id = 1; // 1 as a first bit stands for mirroring of all values
        }else{
            for(int i = 0; i < this->count; i++){
                rank_g = this->translate_up_rank(i)+1;
                if(rank_g >= id_len*32){
                    old_id_len = id_len;
                    id_len = __a_ceil(rank_g/32);
                    hash_id = (unsigned int*)realloc(hash_id, sizeof(unsigned int)*id_len);
                    while(old_id_len < id_len) hash_id[old_id_len++] = 0; // memset with 0
                }
                index = rank_g/32;
                hash_id[index] |= 1 << rank_g % 32;
            }
        }
        return std::pair<unsigned int*,size_t>(hash_id,1);
    }

    group::group(const char* name, int master, group* parent): count(0), members(NULL), members_g(NULL), object_count(0), vacations(NULL), vacant_level(0)
    {
        this->parent = parent;
        this->name = name;
        this->master = master;
        this->parent->children.insert(this);
        group_map(this->name, this);
    }

    void group::add(const int* procs, int count, bool excl){
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

    void group::add(int count, bool excl){
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

    int group::get_vacant(){
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

    void group::add_range(int first, int last, bool excl){
        int count = last - first + 1;
        int procs[count];
        for(int i=0; i<count; i++) procs[i] = first + i;
        this->add(procs, count, excl);
    }

    void group::add_every(int nth, bool excl){
        int count = (int)(this->parent->count / nth);
        int procs[count];
        if(this->parent == NULL) printf("Warning: attempting to modify root group.\n");
        for(int i=0; i < count; i++) procs[i] = (i+1)*nth-1;
        this->add(procs, count, excl);
    }

    void group::add_every(bool(*include)(int k), bool excl){
        int procs[this->parent->count];
        int count = 0;
        for(int i=0; i<this->parent->count; i++){
            if(include(i) == true) procs[count++] = i;
        }
        this->add(procs, count, excl);
    }

    void group::add_every(bool(*include)(int k), bool(*excl)(int k)){
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

    void group::add_intersection(const group* b, int* count){
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

    void group::add_every_intersection(const group* b, int nth, int* count){
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

    void group::add_substraction(const group* b, int* count){
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

    void group::reorder(int* new_ranks){ //  rank i becomes new_ranks[i];
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

    void group::reorder(int(*order)(int r)){
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
    int group::translate_up_rank(int rank, group* parent) const{
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
    int group::translate_down_rank(int rank, group* child) const{
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

    void group::commit(){
        if(this->parent == NULL){ 
            printf("Warning: attempting to commit ambient group.\n");
            return;
        }
        this->id = hash_group_id();
        group** original = mpi_comm_map.get(this->id.first, this->id.second);
       *original == NULL ? *original = this : throw *original;
        MPI_Group_incl(this->parent->mpi_group, this->count, this->members, &this->mpi_group);
        MPI_Comm_create(this->parent->mpi_comm, this->mpi_group, &this->mpi_comm);
        MPI_Group_rank(this->mpi_group, &this->rank);
        if(this->involved()) this->manager = new packet_manager(this);
        this->vacations = (int*)malloc(sizeof(int)*this->count);
        memset(this->vacations, 0, sizeof(int)*this->count);
    }

    int group::get_master(){
        return this->master;
    }

    int group::get_master_g(){
        return translate_up_rank(this->master);
    }

    int group::get_rank(){
        return this->rank;
    }

    void group::spin_loop(){
        assert(this->involved());
        return this->manager->spin_loop();
    }
    void group::spin(){
        assert(this->involved());
        return this->manager->spin();
    }
    packet_manager* group::get_manager(){
        assert(this->involved());
        return this->manager;
    }
    bool group::involved(){
        return (this->rank != MPI_UNDEFINED);
    }
    bool group::is_master(){
        return (this->rank == this->master);
    }
    size_t group::get_size(){
        return this->count;
    }
    const char* group::get_name(){
        return this->name;
    }
    group* group_map(const char* name, group* instance){
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

} }
