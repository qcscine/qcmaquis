#ifndef AMBIENT_GROUPS_GROUP_H
#define AMBIENT_GROUPS_GROUP_H
#include <mpi.h>
#include <set>

#include "ambient/groups/packet_manager.h"
#define ID_TYPE unsigned long long int

namespace ambient{ namespace groups{

    class group
    {
    public:
        group(const char* name, int master, const char* parent);
        group(const char* name, int master, group* parent);
        group(const char* name, int master, MPI_Comm parent); // special constructor for nest group
        void commit();

        void add(const int* procs, int count, bool excl = false);
        void add_range(int first, int last, bool excl = false);
        void add_every(int nth, bool excl = false);
        void add_every(bool(*include)(int k), bool excl = false);
        void add_every(bool(*include)(int k), bool(*excl)(int k));

        void reorder(int* new_ranks);
        void reorder(int(*order)(int r));

        int translate_rank(int rank, group* parent = NULL) const;

        static group* group_map(const char* name, group* instance = NULL);

        const char* name;
        int master;                // master process in this group
        packet_manager* manager;   // group packet manager
        group* parent;             // parent group of processes
        std::set<group*> children;
        MPI_Comm mpi_comm;
        MPI_Group mpi_group;
        int count;                 // number of processes inside group
        int rank;
        ID_TYPE id; // >_<
    private:
        int* members;              // list of member ranks (according to the parent group)
    };

    ID_TYPE get_id();              // generate id of the group

} }
#endif
