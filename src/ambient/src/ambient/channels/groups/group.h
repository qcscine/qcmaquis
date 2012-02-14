#ifndef AMBIENT_GROUPS_GROUP_H
#define AMBIENT_GROUPS_GROUP_H
#include "ambient/utils/singleton.hpp"
#include <stdlib.h>
#include <vector>
#include <set>

namespace ambient { namespace channels {

    class group
    {
    public:
       ~group();
        group(const char* name, int master, group* parent);
        group(const char* name, int master, MPI_Comm parent); // special constructor for nest group
        void commit();

        int get_vacant();   // get the vacant rank for the child creation
        int get_master();   // get rank of the group master
        int get_master_g(); // get translated rank of the group master
        bool occupied();
        bool involved();
        bool is_master();
        void occupy();
        void idle();

        void add(const int* procs, int count, bool excl = false);
        void add(int count, bool excl = false); // loose method
        void add_range(int first, int last, bool excl = false);
        void add_every(int nth, bool excl = false);
        void add_every(bool(*include)(int k), bool excl = false);
        void add_every(bool(*include)(int k), bool(*excl)(int k));
        void add_intersection(const group* b, int* count = NULL);
        void add_every_intersection(const group* b, int nth, int* count = NULL);
        void add_substraction(const group* b, int* count = NULL);
        void reorder(int* new_ranks);
        void reorder(int(*order)(int r));

        int translate_up_rank(int rank, group* parent = NULL) const;
        int translate_down_rank(int rank, group* child) const;

        std::pair<size_t*,size_t> hash_gid();
        int get_rank();
        size_t get_size();
        const char* get_name();
        size_t get_member(size_t i); // returns global children

        int master;                 // master process in this group
        int vacant_level;
        int* vacations;
        group* parent;             // parent group of processes
        std::set<group*> children;
        MPI_Comm mpi_comm;
        MPI_Group mpi_group;
        int count;                 // number of processes inside group
        int rank;
        std::pair<size_t*,size_t> id;
        const char* name;
    private:
        int* members;              // list of member ranks according to the parent group
        int* members_g;            // list of member ranks according to the world group
        bool occupancy;
    };

    group* group_map(const char* name, group* instance = NULL);
    group* group_id_map(std::pair<size_t*,size_t> id, group* instance = NULL);

} }
#endif
