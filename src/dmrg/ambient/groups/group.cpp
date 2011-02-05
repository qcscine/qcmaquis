#include "ambient/ambient.h"
#include "ambient/groups/group.h"
#include "ambient/groups/multirank.h"

namespace ambient{ namespace groups {

    group::group(const char* name, int master, MPI_Comm parent): members(NULL)
    {
        this->parent = NULL;
        this->mpi_comm = parent;
        MPI_Comm_group(this->mpi_comm, &this->mpi_group);
        MPI_Group_size(this->mpi_group, &this->count);
        MPI_Group_rank(this->mpi_group, &this->rank);
        ambient::rank.set( this, this->rank );
        this->name = name;
        this->master = master;
        this->manager = new packet_manager(&this->mpi_comm);
        group_map(name, this);
    }

    group::group(const char* name, int master, group* parent): members(NULL)
    {
        this->parent = parent;
        this->mpi_group = this->parent->mpi_group;
        this->mpi_comm = this->parent->mpi_comm;
        this->name = name;
        this->master = master;
        this->parent->children.insert(this);
        this->manager = new packet_manager(&this->mpi_comm);
        group_map(name, this);
    }

    group::group(const char* name, int master, const char* parent): members(NULL)
    {
        this->parent = group_map(parent);
        this->mpi_group = this->parent->mpi_group;
        this->mpi_comm = this->parent->mpi_comm;
        this->name = name;
        this->master = master;
        this->parent->children.insert(this);
        this->manager = new packet_manager(&this->mpi_comm);
        group_map(name, this);
    }

    void group::add(const int* procs, int count, bool excl){
        this->members = (int*)realloc(this->members, (this->count+count)*sizeof(int));
        memcpy(&(this->members[this->count]), procs, count*sizeof(int));
        if(excl == true){
            // TODO: discard existing groups
        }
        this->count += count;
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

    int group::translate_rank(int rank, group* parent) const{
        int rank_n = rank;
        const group* iterator = this;
        if(!parent) parent = group_map("ambient");
        do{
            if(rank_n >= iterator->count || (iterator->parent == NULL && parent != NULL)){
                printf("Warning: the rank doesn't belongs to this group"); 
                break; 
            }
            rank_n = iterator->members[rank_n];
            iterator = iterator->parent;
        }while(iterator != parent);
        return rank_n;
    }

    void group::commit(){
        if(this->parent == NULL){ 
            printf("Warning: attempting to commit ambient group.\n");
            return;
        }
        MPI_Group_incl(this->parent->mpi_group, this->count, this->members, &this->mpi_group);
        MPI_Comm_create(this->parent->mpi_comm, this->mpi_group, &this->mpi_comm);
        MPI_Group_rank(this->mpi_group, &this->rank);
        ambient::rank.set( this, this->rank );
    }

    group* group::group_map(const char* name, group* instance){
        static std::map<std::string,group*> map;
        if(instance != NULL){ 
            if(map.find(name) != map.end()) printf("Warning: trying to add to groups with the same name\n");
            map.insert(std::pair<const char*,group*>(name,instance));
        }else{
            if(map.find(name) == map.end()) printf("Error: wasn't able to find requested group (%s)\n", name);
            return map.find(name)->second; 
        }
    }

} }
