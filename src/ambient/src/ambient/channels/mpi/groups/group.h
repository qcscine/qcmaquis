#ifndef AMBIENT_CHANNELS_MPI_GROUPS
#define AMBIENT_CHANNELS_MPI_GROUPS

namespace ambient { namespace channels { namespace mpi {

    class group
    {
    public:
        inline ~group();
        inline group(const char* name, int master, group* parent);
        inline group(const char* name, int master, MPI_Comm parent); // special constructor for nest group
        inline void commit();

        inline int get_vacant();   // get the vacant rank for the child creation
        inline int get_master();   // get rank of the group master
        inline int get_master_g(); // get translated rank of the group master
        inline bool occupied();
        inline bool involved();
        inline bool is_master();
        inline void occupy();
        inline void idle();

        inline void add(const int* procs, int count, bool excl = false);
        inline void add(int count, bool excl = false); // loose method
        inline void add_range(int first, int last, bool excl = false);
        inline void add_every(int nth, bool excl = false);
        inline void add_every(bool(*include)(int k), bool excl = false);
        inline void add_every(bool(*include)(int k), bool(*excl)(int k));
        inline void add_intersection(const group* b, int* count = NULL);
        inline void add_every_intersection(const group* b, int nth, int* count = NULL);
        inline void add_substraction(const group* b, int* count = NULL);
        inline void reorder(int* new_ranks);
        inline void reorder(int(*order)(int r));

        inline int translate_up_rank(int rank, group* parent = NULL) const;
        inline int translate_down_rank(int rank, group* child) const;

        inline std::pair<size_t*,size_t> hash_gid();
        inline int get_rank();
        inline size_t get_size();
        inline const char* get_name();
        inline size_t get_member(size_t i); // returns global children

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

    inline group* group_map(const char* name, group* instance = NULL);
    inline group* group_id_map(std::pair<size_t*,size_t> id, group* instance = NULL);

} } }

#endif
