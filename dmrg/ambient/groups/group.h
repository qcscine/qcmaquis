#ifndef AMBIENT_GROUPS_GROUP_H
#define AMBIENT_GROUPS_GROUP_H
#include "ambient/mpi.h"
#include <set>

#include "ambient/groups/packet_manager.h"

namespace ambient{ namespace groups{

    class group;

    class comm_map {
    private: 
        comm_map();                             // constructor is private
        comm_map(comm_map const&);              // copy constructor is private
        comm_map& operator=(comm_map const&);   // assignment operator is private
    public:
        static comm_map& instance();
    public:
        group** get(unsigned int* hash, unsigned int hash_len, int shift = 0) const;
    private:
        mutable std::vector< std::pair<comm_map*, group*> > content;
    };

    class group
    {
    public:
        group(const char* name, int master, group* parent);
        group(const char* name, int master, MPI_Comm parent); // special constructor for nest group
        void commit();

        int vacant_level;
        int* vacations;
        int get_vacant();   // get the vacant rank for the child creation
        int get_master();   // get rank of the group master
        int get_master_g(); // get translated rank of the group master
        bool involved();
        bool is_master();

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

        int master;                // master process in this group
        packet_manager* manager;   // group packet manager
        packet_manager* get_manager();
        void spin();
        void spin_loop();
        group* parent;             // parent group of processes
        std::set<group*> children;
        MPI_Comm mpi_comm;
        MPI_Group mpi_group;
        int count;                 // number of processes inside group
        int rank;
        std::pair<unsigned int*,size_t> id;

        std::pair<unsigned int*,size_t> hash_group_id();
        int get_rank();
        size_t get_size();
        const char* get_name();
    private:
        unsigned int object_count;
        const char* name;
        int* members;              // list of member ranks according to the parent group
        int* members_g;            // list of member ranks according to the world group
    };

    group* group_map(const char* name, group* instance = NULL);

} }
#endif
