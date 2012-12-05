#ifndef AMBIENT_CHANNELS_MPI_GROUPS
#define AMBIENT_CHANNELS_MPI_GROUPS

namespace ambient { namespace channels { namespace mpi {

    class group
    {
    public:
        ~group();
        group(int master, MPI_Comm parent); // global group
        group(int master, group* parent);   // nested group
        void commit();

        size_t get_proc(size_t k);

        void add(const int* procs, int count);
        void add(int count); // loose method
        void add_range(int first, int last);
        void add_every(int nth);
        void add_every(bool(*include)(int k));
        void add_intersection(const group* b, int* count = NULL);
        void add_every_intersection(const group* b, int nth, int* count = NULL);
        void add_substraction(const group* b, int* count = NULL);
        void reorder(int(*permutation)(int r));
        void reorder(int* permutation);

        int get_vacant();    // get the vacant rank for the child creation
        bool occupied();
        void occupy();
        void idle();

        int  vacant_level;
        int* vacations;
        bool occupancy;

        int  id;
        int  rank;
        int  master;                // master process in this group
        int* ranks;                 // member ranks in the parent group
        int* procs;                 // member ranks in the world group
        int  size;                  // number of processes inside group
        int  depth;
        group* parent;              // parent group of processes
        std::set<group*> children;

        MPI_Comm  mpi_comm;
        MPI_Group mpi_group;
    };

    // group* group_map(size_t id, group* instance = NULL);

} } }

#endif
