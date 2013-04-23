#include "ambient/utils/reduce.hpp"

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
        free(this->ranks);
        free(this->procs);
        free(this->vacations);
    }

    inline group::group(int master, MPI_Comm parent)
    : master(master), mpi_comm(parent), id(0), depth(0), parent(NULL),
      ranks(NULL), vacant_level(0), occupancy(false)
    {
        MPI_Comm_group(this->mpi_comm,  &this->mpi_group);
        MPI_Group_size(this->mpi_group, &this->size);
        MPI_Group_rank(this->mpi_group, &this->rank);
        this->procs = (int*)malloc(sizeof(int)*size);
        for(int i = 0; i < size; ++i) this->procs[i] = i;
        this->vacations = (int*)calloc(this->size, sizeof(int));
    }

    inline group::group(int master, group* parent)
    : master(master), parent(parent), depth(parent->depth+1), ranks(NULL), procs(NULL),
      vacations(NULL), vacant_level(0), occupancy(false), size(0)
    {
        this->parent->children.insert(this);
    }

    inline void group::commit(){
        //group_map(id, this);
        MPI_Group_incl(this->parent->mpi_group, this->size, this->ranks, &this->mpi_group);
        MPI_Group_rank(this->mpi_group, &this->rank);
        this->vacations = (int*)calloc(this->size, sizeof(int));
    }

    inline size_t group::get_proc(size_t k){
        assert(k < this->size);
        return this->procs[k];
    }

    // {{{ creation/changing of the group content
    inline void group::add(const int* procs, int count){
        if(count <= 0) return;
        this->ranks = (int*)realloc(this->ranks, (this->size+count)*sizeof(int));
        this->procs = (int*)realloc(this->procs, (this->size+count)*sizeof(int));
        memcpy(&(this->ranks[this->size]), procs, count*sizeof(int));
        for(int i=0; i < count; i++)
            this->procs[this->size+i] = ambient::rank.translate(this->ranks[this->size+i], this->parent);
        this->size += count;
    }

    inline void group::add(int count){
        if(count <= 0) return;
        count = std::min(count, parent->size); // avoiding over-allocation problems
        this->ranks   = (int*)realloc(this->ranks,   (this->size+count)*sizeof(int));
        this->procs = (int*)realloc(this->procs, (this->size+count)*sizeof(int));
        for(int i=0; i < count; i++){ 
            this->ranks  [this->size+i] = this->parent->get_vacant();
            this->procs[this->size+i] = ambient::rank.translate(this->ranks[this->size+i], this->parent);
        }
        this->size += count;
    }

    inline void group::add_range(int first, int last){
        int count = last - first + 1;
        int procs[count];
        for(int i=0; i<count; i++) procs[i] = first + i;
        this->add(procs, count);
    }

    inline void group::add_every(int nth){
        int count = (int)(this->parent->size / nth);
        int procs[count];
        if(this->parent == NULL) printf("Warning: attempting to modify root group.\n");
        for(int i=0; i < count; i++) procs[i] = (i+1)*nth-1;
        this->add(procs, count);
    }

    inline void group::add_every(bool(*include)(int k)){
        int procs[this->parent->size];
        int count = 0;
        for(int i=0; i<this->parent->size; i++)
            if(include(i) == true) procs[count++] = i;
        this->add(procs, count);
    }

    inline void group::add_intersection(const group* b, int* count){
        for(int i=0; i < this->parent->size; i++)
        for(int j=0; j < b->size; j++)
        if(this->parent->procs[i] == b->procs[j]){
            if(count != NULL && (*count)-- <= 0) return;
            this->add(&i, 1);
            break;
        }
    }

    inline void group::add_every_intersection(const group* b, int nth, int* count){
        assert(false);
        int countdown = nth;
        for(int i=0; i < this->parent->size; i++)
        for(int j=0; j < b->size; j++)
        if(this->parent->procs[i] == b->procs[j]){
            if(count != NULL && (*count)-- <= 0) return;
            if(--countdown == 0){ this->add(&i, 1); countdown = nth; }
            break;
        }
    }

    inline void group::add_substraction(const group* b, int* count){
        for(int i=0; i < this->parent->size; i++)
        for(int j=0; j < b->size; j++)
        if(this->parent->procs[i] == b->procs[j]) break;
        else if(j == b->size-1){
            if(count != NULL && (*count)-- <= 0) return;
            this->add(&i, 1);
        }
    }

    inline void group::reorder(int(*permutation)(int r)){
        int* ranks = (int*)malloc(this->size*sizeof(int));
        for(int i=0; i < this->size; i++) ranks[permutation(i)] = this->ranks[i];
        free(this->ranks);
        this->ranks = ranks;

        for(std::set<group*>::iterator it = children.begin(); it != children.end(); ++it)
            for(int i=0; i < (*it)->size; i++)
                (*it)->ranks[i] = permutation((*it)->ranks[i]);
    }

    inline void group::reorder(int* permutation){ //  rank i becomes permutation[i];
        int* ranks = (int*)malloc(this->size*sizeof(int));
        for(int i=0; i < this->size; i++) ranks[permutation[i]] = this->ranks[i];
        free(this->ranks);
        this->ranks = ranks;

        for(std::set<group*>::iterator it = children.begin(); it != children.end(); ++it)
            for(int i=0; i < (*it)->size; i++)
                (*it)->ranks[i] = permutation[(*it)->ranks[i]];
    }
    // }}}

    inline int group::get_vacant(){
        for(int i=0; i < this->size; i++){
            if(this->vacations[i] == this->vacant_level){
                this->vacations[i]++;
                return i;
            }
        }
        this->vacant_level++;
        return this->get_vacant(); 
// note:recursion will be endless in case of 0-group
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

    /*inline group* group_map(size_t id, group* instance){ // timings: 0.36 seconds
        static std::map<size_t,group*> map;
        if(instance != NULL){
            if(map.find(id) != map.end()) map.find(id)->second = instance; // can delete the group here
            else map.insert(std::pair<size_t,group*>(id,instance));
            return instance;
        }
        if(map.find(id) == map.end()) instance = NULL; // wasn't able to find requested group
        else instance = map.find(id)->second;
        return instance;
    }*/

    /*inline std::pair<size_t*,size_t> group::hash_gid(){ // id containing ranks bits
        int index;
        int old_id_len;
        int id_len = 0;
        int rank_g;
        int volume;
        size_t* hash_id = NULL;

        MPI_Comm_size(MPI_COMM_WORLD, &volume);
        if(this->size == volume){
            hash_id = (size_t*)malloc(sizeof(size_t));
            id_len = *hash_id = 1; // 1 as a first bit stands for mirroring of all values
        }else{
            for(int i = 0; i < this->size; i++){
                rank_g = ambient::rank.translate(i, this)+1;
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
    }*/

} } }
